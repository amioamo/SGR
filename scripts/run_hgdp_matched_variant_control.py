"""Within-cohort control on the HGDP harmonized variant set.

Each window uses the sites retained for the HGDP European target. The reference
is the same 1000 Genomes EUR public panel. The only change is the beacon:
the prototype 1000 Genomes cohort, or the HGDP cohort from the cross-dataset run.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from beacon_reconstruction.baselines import pairwise_moments
from beacon_reconstruction.reconstruction import REFERENCE_HAPLOTYPE_RECONSTRUCTION
from beacon_reconstruction.core import ExperimentConfig, BeaconAttackProblem, BeaconEvaluationTruth
from beacon_reconstruction.evaluation import EvaluationResult, evaluate_reconstruction
from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers, matrix_for_samples as _matrix_for_samples
from scripts.lib.sample_splits import load_sample_manifest, load_seed_manifest, split_samples
from scripts.run_hgdp_eur_cross_dataset import (
    METHOD_ORDER,
    METHOD_RUNNERS,
    _format_p,
    _load_hgdp_window,
    _mean,
    _with_config,
    _write_csv,
    _write_json,
)
from scripts.run_template_stage_ablation import participant_similarity_stats
from scripts.lib.paired_statistics import _summarize_deltas


MANIFEST_ID = "1kg_target_on_hgdp_harmonized_sites_v1"
DEFAULT_HGDP_DIR = ROOT / "data" / "external" / "hgdp_1kg"
DEFAULT_MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
DEFAULT_CROSS_METRICS = (
    ROOT / "results" / "block_inheritance_v1" / "hgdp_eur_cross_dataset" / "hgdp_eur_cross_dataset_metrics.csv"
)
DEFAULT_OUT_DIR = ROOT / "results" / "block_inheritance_v1" / "hgdp_matched_variant_control"
MIN_KEPT_VARIANTS = 100
BOOTSTRAP_SEED = 0
BOOTSTRAP_ITERATIONS = 10_000
WITHIN_COHORT = "1kg_target"
CROSS_COHORT = "hgdp_target"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _is_true(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = sorted(
        {(int(row["beacon_size"]), row["method"], row["target_cohort"]) for row in rows},
        key=lambda item: (item[0], METHOD_ORDER.index(item[1]), item[2]),
    )
    out: list[dict[str, Any]] = []
    for beacon_size, method, target_cohort in keys:
        selected = [
            row
            for row in rows
            if int(row["beacon_size"]) == beacon_size
            and row["method"] == method
            and row["target_cohort"] == target_cohort
        ]
        out.append(
            {
                "beacon_size": beacon_size,
                "method": method,
                "target_cohort": target_cohort,
                "n_runs": len(selected),
                "mean_kept_variants": _mean(selected, "kept_variants"),
                "mean_f1": _mean(selected, "f1"),
                "median_f1": float(np.median([float(row["f1"]) for row in selected])),
                "mean_hamming": _mean(selected, "alignment_total_hamming_cost"),
                "mean_similarity": _mean(selected, "mean_similarity"),
                "proportion_similarity_ge_095": _mean(selected, "proportion_similarity_ge_095"),
                "total_source_overlap_count": int(sum(int(row["source_overlap_count"]) for row in selected)),
                "invalid_run_count": int(sum(_is_true(row["invalid_run"]) for row in selected)),
            }
        )
    return out


def _window_means(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = sorted(
        {
            (row["window_id"], int(row["beacon_size"]), row["method"], int(row["kept_variants"]))
            for row in rows
            if row["target_cohort"] == WITHIN_COHORT
        }
    )
    out: list[dict[str, Any]] = []
    for window_id, beacon_size, method, kept_variants in keys:
        within = [
            row
            for row in rows
            if row["window_id"] == window_id
            and int(row["beacon_size"]) == beacon_size
            and row["method"] == method
            and row["target_cohort"] == WITHIN_COHORT
        ]
        cross = [
            row
            for row in rows
            if row["window_id"] == window_id
            and int(row["beacon_size"]) == beacon_size
            and row["method"] == method
            and row["target_cohort"] == CROSS_COHORT
        ]
        out.append(
            {
                "window_id": window_id,
                "beacon_size": beacon_size,
                "method": method,
                "kept_variants": kept_variants,
                "n_replicates": len(within),
                "within_f1": _mean(within, "f1"),
                "cross_f1": _mean(cross, "f1"),
                "f1_cross_minus_within": _mean(cross, "f1") - _mean(within, "f1"),
            }
        )
    return out


def _paired_inference(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {
        (int(row["beacon_size"]), row["window_id"], str(row["replicate_id"]), row["method"], row["target_cohort"]): row
        for row in rows
    }
    units = sorted(
        {
            (int(row["beacon_size"]), row["window_id"], str(row["replicate_id"]), row["method"])
            for row in rows
            if row["target_cohort"] == WITHIN_COHORT
        }
    )
    grouped: dict[tuple[int, str], list[tuple[str, str]]] = {}
    for beacon_size, window_id, replicate_id, method in units:
        grouped.setdefault((beacon_size, method), []).append((window_id, replicate_id))
    inference: list[dict[str, Any]] = []
    for (beacon_size, method), pairs in sorted(grouped.items(), key=lambda item: (item[0][0], METHOD_ORDER.index(item[0][1]))):
        for metric in ("f1", "proportion_similarity_ge_095"):
            deltas: list[float] = []
            baseline_values: list[float] = []
            method_values: list[float] = []
            for window_id, replicate_id in pairs:
                within = by_key[(beacon_size, window_id, replicate_id, method, WITHIN_COHORT)]
                cross = by_key[(beacon_size, window_id, replicate_id, method, CROSS_COHORT)]
                deltas.append(float(cross[metric]) - float(within[metric]))
                baseline_values.append(float(within[metric]))
                method_values.append(float(cross[metric]))
            summary = _summarize_deltas(
                comparison="hgdp_target_minus_1kg_target",
                metric=metric,
                deltas=deltas,
                baseline_values=baseline_values,
                method_values=method_values,
                higher_is_better=True,
                bootstrap_seed=BOOTSTRAP_SEED + beacon_size + (0 if method == REFERENCE_HAPLOTYPE_RECONSTRUCTION else 17),
                bootstrap_iterations=BOOTSTRAP_ITERATIONS,
            )
            summary["beacon_size"] = beacon_size
            summary["method"] = method
            inference.append(summary)
    return inference


def _write_report(
    path: Path,
    *,
    aggregate_rows: list[dict[str, Any]],
    window_rows: list[dict[str, Any]],
    inference_rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> None:
    lines = [
        "# Matched-Variant Within-Cohort Control",
        "",
        "The variant set in each window is the set retained for the HGDP European target.",
        "Both arms use the same 1000 Genomes EUR public reference. The within-cohort arm",
        "replaces the HGDP beacon with the prototype 1000 Genomes beacon.",
        "",
        "## Configuration",
        "",
        f"- beacon sizes: {', '.join(str(value) for value in config['beacon_sizes'])}",
        f"- reference panel size: {config['panel_size']}",
        f"- replicates: {config['max_replicates']}",
        f"- minimum kept variants: {config['min_kept_variants']}",
        "",
        "## Carrier F1 on the Shared Sites",
        "",
        "| N | method | target | n | mean F1 | median F1 | mean Hamming | Sim ≥ 0.95 | invalid |",
        "|---:|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['beacon_size']} | {row['method']} | {row['target_cohort']} | {row['n_runs']} | "
            f"{row['mean_f1']:.6f} | {row['median_f1']:.6f} | {row['mean_hamming']:.2f} | "
            f"{row['proportion_similarity_ge_095']:.6f} | {row['invalid_run_count']} |"
        )
    lines.extend(
        [
            "",
            "## Paired Difference, HGDP Target minus 1KG Target",
            "",
            "Positive means the HGDP beacon is reconstructed more accurately than the 1KG beacon",
            "on the same sites and the same reference panel.",
            "",
            "| N | method | metric | n | mean within | mean cross | mean delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |",
            "|---:|---|---|---:|---:|---:|---:|---|---:|---:|---:|",
        ]
    )
    for row in inference_rows:
        ci = f"[{row['ci95_low']:.6f}, {row['ci95_high']:.6f}]"
        wtl = f"{row['wins']}/{row['ties']}/{row['losses']}"
        lines.append(
            f"| {row['beacon_size']} | {row['method']} | {row['metric']} | {row['n_pairs']} | "
            f"{row['mean_baseline']:.6f} | {row['mean_method']:.6f} | {row['mean_delta']:.6f} | {ci} | "
            f"{wtl} | {_format_p(float(row['sign_test_p']))} | {_format_p(float(row['wilcoxon_p']))} |"
        )
    lines.extend(
        [
            "",
            "## Full Method by Window",
            "",
            "| window | sites | N | within F1 | cross F1 | cross − within |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in window_rows:
        if row["method"] != REFERENCE_HAPLOTYPE_RECONSTRUCTION:
            continue
        lines.append(
            f"| {row['window_id']} | {row['kept_variants']} | {row['beacon_size']} | "
            f"{row['within_f1']:.6f} | {row['cross_f1']:.6f} | {row['f1_cross_minus_within']:.6f} |"
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "This comparison fixes the harmonized site set. It does not restore the original",
            "500-SNP windows. The two beacons are different people, so the paired difference",
            "mixes cohort shift with the particular beacon draw.",
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def _cross_rows(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in _read_csv(path):
        if int(row["kept_variants"]) < MIN_KEPT_VARIANTS:
            continue
        rows.append(
            {
                "window_id": row["window_id"],
                "replicate_id": row["replicate_id"],
                "beacon_size": int(row["beacon_size"]),
                "method": row["method"],
                "target_cohort": CROSS_COHORT,
                "reference_cohort": "1kg_eur",
                "kept_variants": int(row["kept_variants"]),
                "precision": float(row["precision"]),
                "recall": float(row["recall"]),
                "f1": float(row["f1"]),
                "alignment_total_hamming_cost": float(row["alignment_total_hamming_cost"]),
                "count_mae": float(row["count_mae"]),
                "max_count_error": float(row["max_count_error"]),
                "mean_similarity": float(row["mean_similarity"]),
                "proportion_similarity_ge_095": float(row["proportion_similarity_ge_095"]),
                "source_overlap_count": int(row["source_overlap_count"]),
                "invalid_run": _is_true(row["invalid_run"]),
                "invalid_reasons": row["invalid_reasons"],
            }
        )
    if not rows:
        raise ValueError(f"No cross-cohort rows in {path}.")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hgdp-dir", type=Path, default=DEFAULT_HGDP_DIR)
    parser.add_argument("--multiwindow-dir", type=Path, default=DEFAULT_MULTIWINDOW_DIR)
    parser.add_argument("--cross-metrics", type=Path, default=DEFAULT_CROSS_METRICS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--beacon-sizes", default="25,50,100")
    args = parser.parse_args()

    windows = _read_csv(args.multiwindow_dir / "multiwindow_manifest.csv")[: max(0, args.max_windows)]
    overlap_by_window = {row["window_id"]: row for row in _read_csv(args.hgdp_dir / "hgdp_eur_window_overlap.csv")}
    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]
    beacon_sizes = [int(value) for value in args.beacon_sizes.split(",") if value.strip()]
    rows = _cross_rows(args.cross_metrics)
    for window in windows:
        window_id = window["window_id"]
        overlap = overlap_by_window[window_id]
        if int(overlap["kept_variants"]) < MIN_KEPT_VARIANTS:
            continue
        _hgdp_ids, hgdp_variants, hgdp_positions, _hgdp_matrix = _load_hgdp_window(
            args.hgdp_dir / "windows" / f"{window_id}.npz"
        )
        if len(hgdp_variants) != int(overlap["kept_variants"]):
            raise ValueError(f"{window_id} npz has {len(hgdp_variants)} variants, overlap says {overlap['kept_variants']}.")
        kg_ids, kg_variants, _positions, _alleles, kg_matrix = _load_raw_as_target_carriers(Path(window["raw_path"]))
        kg_index = {variant_id: index for index, variant_id in enumerate(kg_variants)}
        missing_variants = [variant_id for variant_id in hgdp_variants if variant_id not in kg_index]
        if missing_variants:
            raise ValueError(f"{window_id} is missing {len(missing_variants)} harmonized variants in the 1KG matrix.")
        kg_matrix = kg_matrix[:, [kg_index[variant_id] for variant_id in hgdp_variants]]
        for beacon_size in beacon_sizes:
            for replicate_id in replicate_ids:
                seeds = seeds_by_replicate[replicate_id]
                split = split_samples(sample_rows, replicate_id, beacon_size, shared_reference_panel=True)
                reference_ids = split["public_ref"]
                beacon_ids = split["beacon"]
                overlap_count = len(set(beacon_ids).intersection(reference_ids))
                if overlap_count:
                    raise ValueError(f"Beacon/reference overlap at {window_id} replicate {replicate_id}.")
                beacon_matrix = _matrix_for_samples(kg_ids, kg_matrix, beacon_ids)
                reference_matrix = _matrix_for_samples(kg_ids, kg_matrix, reference_ids)
                target_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                reference_pairwise = pairwise_moments(reference_matrix)
                split_id = f"{window_id}_1kg_matched_sites_replicate_{replicate_id}_N{beacon_size}"
                problem = BeaconAttackProblem(
                    manifest_id=MANIFEST_ID,
                    split_id=split_id,
                    beacon_size=beacon_size,
                    variant_ids=hgdp_variants,
                    positions=hgdp_positions,
                    query_responses=target_counts,
                    target_counts=target_counts,
                    public_ld=reference_pairwise,
                    config=ExperimentConfig(
                        initialization_seed=seeds["initialization_seed"],
                        optimizer_seed=seeds["optimizer_seed"],
                        group_selection_seed=seeds["group_selection_seed"],
                        metadata={"window_id": window_id, "experiment": "hgdp_matched_variant_control"},
                    ),
                )
                truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(beacon_ids))
                runtime_config = {
                    "window_id": window_id,
                    "replicate_id": replicate_id,
                    "beacon_size": beacon_size,
                    "kept_variants": len(hgdp_variants),
                    "target_cohort": WITHIN_COHORT,
                    "reference_cohort": "1kg_eur",
                }
                for method in METHOD_ORDER:
                    result = _with_config(
                        METHOD_RUNNERS[method](problem, reference_matrix, metadata=runtime_config),
                        config=runtime_config,
                        optimizer_seed=seeds["optimizer_seed"],
                    )
                    evaluation: EvaluationResult = evaluate_reconstruction(
                        result,
                        truth,
                        target_counts=target_counts,
                        reference_pairwise=reference_pairwise,
                        source_sets={"beacon": beacon_ids, "template_reference": reference_ids},
                        accepted_manifest_ids=None,
                    )
                    similarity = participant_similarity_stats(beacon_matrix, result.reconstructed_binary)
                    metrics = evaluation.metrics
                    rows.append(
                        {
                            "window_id": window_id,
                            "replicate_id": replicate_id,
                            "beacon_size": beacon_size,
                            "method": method,
                            "target_cohort": WITHIN_COHORT,
                            "reference_cohort": "1kg_eur",
                            "kept_variants": len(hgdp_variants),
                            "precision": metrics.precision,
                            "recall": metrics.recall,
                            "f1": metrics.f1,
                            "alignment_total_hamming_cost": metrics.alignment_total_hamming_cost,
                            "count_mae": metrics.count_mae,
                            "max_count_error": metrics.max_count_error,
                            "mean_similarity": similarity["mean_similarity"],
                            "proportion_similarity_ge_095": similarity["proportion_similarity_ge_095"],
                            "source_overlap_count": overlap_count,
                            "invalid_run": evaluation.invalid_run.is_invalid,
                            "invalid_reasons": ";".join(evaluation.invalid_run.reasons),
                        }
                    )
                print(f"finished {split_id}", flush=True)

    aggregate_rows = _aggregate(rows)
    window_rows = _window_means(rows)
    inference_rows = _paired_inference(rows)
    config = {
        "beacon_sizes": beacon_sizes,
        "panel_size": 100,
        "max_replicates": args.max_replicates,
        "min_kept_variants": MIN_KEPT_VARIANTS,
        "cross_metrics": str(args.cross_metrics),
        "site_rule": "hgdp_harmonized_sites_called_in_every_retained_sample",
    }
    out_dir = args.out_dir
    _write_csv(out_dir / "hgdp_matched_variant_control_metrics.csv", rows)
    _write_csv(out_dir / "hgdp_matched_variant_control_aggregate.csv", aggregate_rows)
    _write_csv(out_dir / "hgdp_matched_variant_control_by_window.csv", window_rows)
    _write_csv(out_dir / "hgdp_matched_variant_control_inference.csv", inference_rows)
    _write_json(out_dir / "hgdp_matched_variant_control_manifest.json", config)
    _write_report(
        out_dir / "hgdp_matched_variant_control_report.md",
        aggregate_rows=aggregate_rows,
        window_rows=window_rows,
        inference_rows=inference_rows,
        config=config,
    )


if __name__ == "__main__":
    main()
