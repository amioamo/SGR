"""Reconstruct an HGDP European beacon with a 1000 Genomes EUR reference.

The beacon is drawn from the HGDP European subset. The reference is the
paper's matched 1000 Genomes EUR public panel for the same replicate. The two
cohorts come from different studies.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
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
from beacon_reconstruction.reconstruction import REFERENCE_ONLY, REFERENCE_HAPLOTYPE_RECONSTRUCTION, reference_only, reference_haplotype_reconstruction
from beacon_reconstruction.core import ExperimentConfig, BeaconAttackProblem, BeaconEvaluationTruth, ReconstructionResult
from beacon_reconstruction.evaluation import EvaluationResult, evaluate_reconstruction
from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers, matrix_for_samples as _matrix_for_samples
from scripts.lib.sample_splits import load_sample_manifest, load_seed_manifest, split_samples
from scripts.run_template_stage_ablation import participant_similarity_stats
from scripts.lib.paired_statistics import _summarize_deltas


MANIFEST_ID = "hgdp_eur_target_1kg_eur_reference_v1"
DEFAULT_HGDP_DIR = ROOT / "data" / "external" / "hgdp_1kg"
DEFAULT_MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
DEFAULT_OUT_DIR = ROOT / "results" / "block_inheritance_v1" / "hgdp_eur_cross_dataset"
METHOD_ORDER = (REFERENCE_ONLY, REFERENCE_HAPLOTYPE_RECONSTRUCTION)
METHOD_RUNNERS = {
    REFERENCE_ONLY: reference_only,
    REFERENCE_HAPLOTYPE_RECONSTRUCTION: reference_haplotype_reconstruction,
}
BOOTSTRAP_SEED = 0
BOOTSTRAP_ITERATIONS = 10_000


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _serialize_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, dict):
        return {str(key): _serialize_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize_value(item) for item in value]
    return value


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}.")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_serialize_value(payload), indent=2, sort_keys=True) + "\n")


def _mean(rows: list[dict[str, Any]], key: str) -> float:
    values = np.asarray([float(row[key]) for row in rows], dtype=float)
    return float(np.nanmean(values)) if values.size else float("nan")


def _format_p(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    if value < 1e-4:
        return f"{value:.2e}"
    return f"{value:.4f}"


def _load_hgdp_window(path: Path) -> tuple[list[str], tuple[str, ...], tuple[int, ...], np.ndarray]:
    with np.load(path, allow_pickle=False) as payload:
        sample_ids = [str(sample_id) for sample_id in payload["sample_ids"].tolist()]
        variant_ids = tuple(str(variant_id) for variant_id in payload["variant_ids"].tolist())
        positions = tuple(int(position) for position in payload["positions"].tolist())
        carriers = payload["carriers"].astype(np.int8, copy=True)
    if carriers.shape != (len(sample_ids), len(variant_ids)):
        raise ValueError(f"Unexpected HGDP matrix shape in {path}: {carriers.shape}.")
    return sample_ids, variant_ids, positions, carriers


def _nested_beacon_ids(sample_ids: list[str], beacon_size: int, seed: int) -> list[str]:
    if beacon_size > len(sample_ids):
        raise ValueError(f"Beacon size {beacon_size} exceeds the HGDP cohort of {len(sample_ids)}.")
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(sample_ids))
    return [sample_ids[int(index)] for index in order[:beacon_size]]


def _with_config(
    result: ReconstructionResult,
    *,
    config: dict[str, Any],
    optimizer_seed: int,
) -> ReconstructionResult:
    seeds = dict(result.seeds)
    seeds.setdefault("optimizer_seed", int(optimizer_seed))
    metadata = dict(result.metadata)
    metadata["config"] = config
    return ReconstructionResult(
        manifest_id=result.manifest_id,
        split_id=result.split_id,
        prior_id=result.prior_id,
        method_name=result.method_name,
        seeds=seeds,
        reconstructed_binary=result.reconstructed_binary,
        metrics=result.metrics,
        metadata=metadata,
    )


def _aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = sorted(
        {(int(row["beacon_size"]), row["method"]) for row in rows},
        key=lambda item: (item[0], METHOD_ORDER.index(item[1])),
    )
    out: list[dict[str, Any]] = []
    for beacon_size, method in keys:
        selected = [row for row in rows if int(row["beacon_size"]) == beacon_size and row["method"] == method]
        out.append(
            {
                "beacon_size": beacon_size,
                "method": method,
                "n_runs": len(selected),
                "mean_f1": _mean(selected, "f1"),
                "median_f1": float(np.median([float(row["f1"]) for row in selected])),
                "mean_kept_variants": _mean(selected, "kept_variants"),
                "mean_hamming": _mean(selected, "alignment_total_hamming_cost"),
                "mean_count_mae": _mean(selected, "count_mae"),
                "mean_max_count_error": _mean(selected, "max_count_error"),
                "mean_similarity": _mean(selected, "mean_similarity"),
                "proportion_similarity_ge_095": _mean(selected, "proportion_similarity_ge_095"),
                "total_source_overlap_count": int(sum(int(row["source_overlap_count"]) for row in selected)),
                "invalid_run_count": int(sum(str(row["invalid_run"]).lower() in {"1", "true"} for row in selected)),
            }
        )
    return out


def _paired_inference(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {
        (int(row["beacon_size"]), row["window_id"], row["replicate_id"], row["method"]): row for row in rows
    }
    units = sorted({(int(row["beacon_size"]), row["window_id"], row["replicate_id"]) for row in rows})
    grouped: dict[int, list[tuple[str, str]]] = {}
    for beacon_size, window_id, replicate_id in units:
        grouped.setdefault(beacon_size, []).append((window_id, replicate_id))
    inference: list[dict[str, Any]] = []
    for beacon_size, pairs in sorted(grouped.items()):
        for metric in ("f1", "proportion_similarity_ge_095"):
            deltas: list[float] = []
            baseline_values: list[float] = []
            method_values: list[float] = []
            for window_id, replicate_id in pairs:
                baseline = by_key[(beacon_size, window_id, replicate_id, REFERENCE_ONLY)]
                method = by_key[(beacon_size, window_id, replicate_id, REFERENCE_HAPLOTYPE_RECONSTRUCTION)]
                deltas.append(float(method[metric]) - float(baseline[metric]))
                baseline_values.append(float(baseline[metric]))
                method_values.append(float(method[metric]))
            summary = _summarize_deltas(
                comparison="reference_haplotype_reconstruction_minus_reference_only",
                metric=metric,
                deltas=deltas,
                baseline_values=baseline_values,
                method_values=method_values,
                higher_is_better=True,
                bootstrap_seed=BOOTSTRAP_SEED + beacon_size,
                bootstrap_iterations=BOOTSTRAP_ITERATIONS,
            )
            summary["beacon_size"] = beacon_size
            inference.append(summary)
    return inference


def _write_report(
    path: Path,
    *,
    aggregate_rows: list[dict[str, Any]],
    inference_rows: list[dict[str, Any]],
    overlap_rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> None:
    lines = [
        "# HGDP European Target, 1000 Genomes EUR Reference",
        "",
        "The beacon is a sample of HGDP Europeans from the gnomAD v3.1.2 HGDP+1KG",
        "callset. The reference is the matched 1000 Genomes EUR public panel used in",
        "the main experiment. Sites are restricted to variants present in both, with the",
        "1000 Genomes target allele, and called in every retained HGDP sample.",
        "Samples below a 95% call rate on the matched sites are removed before that filter.",
        "",
        "## Configuration",
        "",
        f"- HGDP European samples: {config['hgdp_sample_count']}",
        f"- beacon sizes: {', '.join(str(value) for value in config['beacon_sizes'])}",
        f"- reference panel size: {config['panel_size']}",
        f"- windows requested: {config['max_windows']}",
        f"- windows analyzed: {sum(int(row['kept_variants']) >= 100 for row in overlap_rows)}",
        f"- replicates: {config['max_replicates']}",
        "",
        "## Variant Overlap",
        "",
        "| window | 1KG variants | kept | exact | flipped | missing site | dropped missing GT |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in overlap_rows:
        lines.append(
            f"| {row['window_id']} | {row['target_variants']} | {row['kept_variants']} | {row['exact']} | "
            f"{row['flipped']} | {row['missing_site']} | {row['dropped_missing_genotype']} |"
        )
    lines.extend(
        [
            "",
            "## Carrier F1",
            "",
            "| N | method | n | mean F1 | median F1 | mean Hamming | mean Sim | Sim ≥ 0.95 | invalid |",
            "|---:|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in aggregate_rows:
        lines.append(
            f"| {row['beacon_size']} | {row['method']} | {row['n_runs']} | {row['mean_f1']:.6f} | "
            f"{row['median_f1']:.6f} | {row['mean_hamming']:.2f} | {row['mean_similarity']:.6f} | "
            f"{row['proportion_similarity_ge_095']:.6f} | {row['invalid_run_count']} |"
        )
    lines.extend(
        [
            "",
            "## Paired Difference, Full Method minus Reference-only",
            "",
            "| N | metric | n | mean delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |",
            "|---:|---|---:|---:|---|---:|---:|---:|",
        ]
    )
    for row in inference_rows:
        ci = f"[{row['ci95_low']:.6f}, {row['ci95_high']:.6f}]"
        wtl = f"{row['wins']}/{row['ties']}/{row['losses']}"
        lines.append(
            f"| {row['beacon_size']} | {row['metric']} | {row['n_pairs']} | {row['mean_delta']:.6f} | {ci} | "
            f"{wtl} | {_format_p(float(row['sign_test_p']))} | {_format_p(float(row['wilcoxon_p']))} |"
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "HGDP and 1000 Genomes are different sample collections. N=250 is not run because",
            "the retained HGDP European subset has 146 individuals. Windows with fewer than",
            "100 harmonized variants are omitted; chr20:29757801-29829649 has no variants in",
            "this gnomAD VCF. Reference-only draws rows from the 1000 Genomes panel and does",
            "not use the beacon counts.",
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hgdp-dir", type=Path, default=DEFAULT_HGDP_DIR)
    parser.add_argument("--multiwindow-dir", type=Path, default=DEFAULT_MULTIWINDOW_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--beacon-sizes", default="25,50,100")
    args = parser.parse_args()

    windows = _read_csv(args.multiwindow_dir / "multiwindow_manifest.csv")[: max(0, args.max_windows)]
    overlap_rows = _read_csv(args.hgdp_dir / "hgdp_eur_window_overlap.csv")
    overlap_by_window = {row["window_id"]: row for row in overlap_rows}
    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]
    beacon_sizes = [int(value) for value in args.beacon_sizes.split(",") if value.strip()]
    hgdp_samples = [
        row
        for row in _read_csv(args.hgdp_dir / "hgdp_eur_samples.csv")
        if str(row.get("retained", "true")).lower() == "true"
    ]
    config = {
        "hgdp_sample_count": len(hgdp_samples),
        "beacon_sizes": beacon_sizes,
        "panel_size": 100,
        "max_windows": args.max_windows,
        "max_replicates": args.max_replicates,
        "target_cohort": "HGDP_EUR",
        "reference_cohort": "1000G_EUR_public_panel",
        "hgdp_dir": str(args.hgdp_dir),
        "out_dir": str(args.out_dir),
    }
    rows: list[dict[str, Any]] = []
    for window in windows:
        window_id = window["window_id"]
        overlap = overlap_by_window[window_id]
        if int(overlap["kept_variants"]) < 100:
            print(f"skip {window_id}: {overlap['kept_variants']} harmonized variants", flush=True)
            continue
        hgdp_ids, hgdp_variants, hgdp_positions, hgdp_matrix = _load_hgdp_window(
            args.hgdp_dir / "windows" / f"{window_id}.npz"
        )
        kg_ids, kg_variants, _kg_positions, _kg_alleles, kg_matrix = _load_raw_as_target_carriers(
            Path(window["raw_path"])
        )
        kg_index = {variant_id: index for index, variant_id in enumerate(kg_variants)}
        if any(variant_id not in kg_index for variant_id in hgdp_variants):
            raise ValueError(f"{window_id} has an HGDP variant absent from the 1000 Genomes window.")
        kg_matrix = kg_matrix[:, [kg_index[variant_id] for variant_id in hgdp_variants]]
        shared_ids = set(hgdp_ids).intersection(kg_ids)
        if shared_ids:
            raise ValueError(f"{window_id} has {len(shared_ids)} sample IDs in both cohorts.")
        for beacon_size in beacon_sizes:
            for replicate_id in replicate_ids:
                seeds = seeds_by_replicate[replicate_id]
                split = split_samples(sample_rows, replicate_id, beacon_size, shared_reference_panel=True)
                reference_ids = split["public_ref"]
                if len(reference_ids) != 100:
                    raise ValueError(f"Expected 100 reference samples, found {len(reference_ids)}.")
                beacon_ids = _nested_beacon_ids(hgdp_ids, beacon_size, seeds["cohort_seed"])
                beacon_matrix = _matrix_for_samples(hgdp_ids, hgdp_matrix, beacon_ids)
                reference_matrix = _matrix_for_samples(kg_ids, kg_matrix, reference_ids)
                target_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                reference_pairwise = pairwise_moments(reference_matrix)
                split_id = f"{window_id}_hgdp_eur_replicate_{replicate_id}_N{beacon_size}"
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
                        metadata={"window_id": window_id, "experiment": "hgdp_eur_cross_dataset"},
                    ),
                )
                truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(beacon_ids))
                runtime_config = {
                    "window_id": window_id,
                    "replicate_id": replicate_id,
                    "beacon_size": beacon_size,
                    "kept_variants": len(hgdp_variants),
                    "target_cohort": "HGDP_EUR",
                    "reference_cohort": "1000G_EUR",
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
                            "kept_variants": len(hgdp_variants),
                            "hgdp_sample_count": len(hgdp_ids),
                            "panel_size": len(reference_ids),
                            "precision": metrics.precision,
                            "recall": metrics.recall,
                            "f1": metrics.f1,
                            "alignment_total_hamming_cost": metrics.alignment_total_hamming_cost,
                            "count_mae": metrics.count_mae,
                            "max_count_error": metrics.max_count_error,
                            "mean_similarity": similarity["mean_similarity"],
                            "proportion_similarity_ge_095": similarity["proportion_similarity_ge_095"],
                            "source_overlap_count": 0,
                            "invalid_run": evaluation.invalid_run.is_invalid,
                            "invalid_reasons": ";".join(evaluation.invalid_run.reasons),
                        }
                    )
                print(f"finished {split_id}", flush=True)

    if not rows:
        raise ValueError("No HGDP windows had enough harmonized variants to run.")
    aggregate_rows = _aggregate(rows)
    inference_rows = _paired_inference(rows)
    out_dir = args.out_dir
    _write_csv(out_dir / "hgdp_eur_cross_dataset_metrics.csv", rows)
    _write_csv(out_dir / "hgdp_eur_cross_dataset_aggregate.csv", aggregate_rows)
    _write_csv(out_dir / "hgdp_eur_cross_dataset_inference.csv", inference_rows)
    _write_json(out_dir / "hgdp_eur_cross_dataset_manifest.json", config)
    _write_report(
        out_dir / "hgdp_eur_cross_dataset_report.md",
        aggregate_rows=aggregate_rows,
        inference_rows=inference_rows,
        overlap_rows=overlap_rows,
        config=config,
    )


if __name__ == "__main__":
    main()
