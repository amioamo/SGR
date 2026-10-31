"""Run Count, LD/AF optimization, and SGR on the HGDP European beacon.

The beacon, sites, and 1000 Genomes EUR reference match the cross-dataset run.
Count is independent count assignment. LD/AF is the Saleem alternating optimizer
with its published oracle initialization. SGR is reference haplotype reconstruction.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from beacon_reconstruction.baselines import INDEPENDENT_COUNT, independent_count
from beacon_reconstruction.baselines import pairwise_moments
from beacon_reconstruction.reconstruction import REFERENCE_HAPLOTYPE_RECONSTRUCTION, reference_haplotype_reconstruction
from beacon_reconstruction.core import ExperimentConfig, BeaconAttackProblem, BeaconEvaluationTruth, ReconstructionResult
from beacon_reconstruction.evaluation import EvaluationResult, evaluate_reconstruction
from beacon_reconstruction.ld_af_baseline import SALEEM_RECONSTRUCTION_ATTACK, saleem_reconstruction_attack
from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers, matrix_for_samples as _matrix_for_samples
from scripts.lib.sample_splits import load_sample_manifest, load_seed_manifest, split_samples
from scripts.run_hgdp_eur_cross_dataset import (
    _format_p,
    _load_hgdp_window,
    _mean,
    _nested_beacon_ids,
    _with_config,
    _write_csv,
    _write_json,
)
from scripts.run_template_stage_ablation import participant_similarity_stats
from scripts.lib.paired_statistics import _summarize_deltas


MANIFEST_ID = "hgdp_eur_target_1kg_eur_reference_v1"
DEFAULT_HGDP_DIR = ROOT / "data" / "external" / "hgdp_1kg"
DEFAULT_MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
DEFAULT_SGR_METRICS = (
    ROOT / "results" / "block_inheritance_v1" / "hgdp_eur_cross_dataset" / "hgdp_eur_cross_dataset_metrics.csv"
)
DEFAULT_OUT_DIR = ROOT / "results" / "block_inheritance_v1" / "hgdp_eur_paper_methods"
METHOD_ORDER = (INDEPENDENT_COUNT, SALEEM_RECONSTRUCTION_ATTACK, REFERENCE_HAPLOTYPE_RECONSTRUCTION)
METHOD_LABELS = {
    INDEPENDENT_COUNT: "Count-based",
    SALEEM_RECONSTRUCTION_ATTACK: "LD/AF optimization",
    REFERENCE_HAPLOTYPE_RECONSTRUCTION: "SGR",
}
BOOTSTRAP_SEED = 0
BOOTSTRAP_ITERATIONS = 10_000
MIN_KEPT_VARIANTS = 100


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _run_method(
    method: str,
    problem: BeaconAttackProblem,
    reference_matrix: np.ndarray,
    beacon_matrix: np.ndarray,
    reference_pairwise: np.ndarray,
    runtime_config: dict[str, Any],
) -> ReconstructionResult:
    if method == INDEPENDENT_COUNT:
        return independent_count(problem, reference_pairwise=reference_pairwise)
    if method == SALEEM_RECONSTRUCTION_ATTACK:
        return saleem_reconstruction_attack(
            problem,
            reference_matrix,
            ground_truth_beacon=beacon_matrix,
            corr_epoch=1000,
            freq_epoch=500,
            max_iterations=3000,
            flip_threshold=10,
            min_outer_iterations=1,
            learning_rate=0.001,
            device="cpu",
        )
    if method == REFERENCE_HAPLOTYPE_RECONSTRUCTION:
        return reference_haplotype_reconstruction(problem, reference_matrix, metadata=runtime_config)
    raise ValueError(f"Unsupported method {method}.")


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
                "method_label": METHOD_LABELS[method],
                "n_runs": len(selected),
                "mean_f1": _mean(selected, "f1"),
                "sd_f1": float(np.std([float(row["f1"]) for row in selected], ddof=1)),
                "median_f1": float(np.median([float(row["f1"]) for row in selected])),
                "mean_hamming": _mean(selected, "alignment_total_hamming_cost"),
                "mean_similarity": _mean(selected, "mean_similarity"),
                "proportion_similarity_ge_095": _mean(selected, "proportion_similarity_ge_095"),
                "mean_wall_time_sec": _mean(selected, "wall_time_sec"),
                "invalid_run_count": int(sum(str(row["invalid_run"]).lower() in {"1", "true"} for row in selected)),
            }
        )
    return out


def _paired_inference(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {
        (int(row["beacon_size"]), row["window_id"], str(row["replicate_id"]), row["method"]): row for row in rows
    }
    units = sorted({(int(row["beacon_size"]), row["window_id"], str(row["replicate_id"])) for row in rows})
    grouped: dict[int, list[tuple[str, str]]] = {}
    for beacon_size, window_id, replicate_id in units:
        grouped.setdefault(beacon_size, []).append((window_id, replicate_id))
    inference: list[dict[str, Any]] = []
    for beacon_size, pairs in sorted(grouped.items()):
        for baseline in (INDEPENDENT_COUNT, SALEEM_RECONSTRUCTION_ATTACK):
            for metric in ("f1", "proportion_similarity_ge_095"):
                deltas: list[float] = []
                baseline_values: list[float] = []
                method_values: list[float] = []
                for window_id, replicate_id in pairs:
                    base = by_key[(beacon_size, window_id, replicate_id, baseline)]
                    method = by_key[(beacon_size, window_id, replicate_id, REFERENCE_HAPLOTYPE_RECONSTRUCTION)]
                    deltas.append(float(method[metric]) - float(base[metric]))
                    baseline_values.append(float(base[metric]))
                    method_values.append(float(method[metric]))
                summary = _summarize_deltas(
                    comparison=f"sgr_minus_{baseline}",
                    metric=metric,
                    deltas=deltas,
                    baseline_values=baseline_values,
                    method_values=method_values,
                    higher_is_better=True,
                    bootstrap_seed=BOOTSTRAP_SEED + beacon_size + (0 if baseline == INDEPENDENT_COUNT else 19),
                    bootstrap_iterations=BOOTSTRAP_ITERATIONS,
                )
                summary["beacon_size"] = beacon_size
                summary["baseline"] = baseline
                inference.append(summary)
    return inference


def _write_report(
    path: Path,
    *,
    aggregate_rows: list[dict[str, Any]],
    inference_rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> None:
    lines = [
        "# HGDP European Target: Count, LD/AF, and SGR",
        "",
        "The beacon is the retained HGDP European cohort. The reference is the 1000 Genomes",
        "EUR public panel. Sites are the harmonized variants used in the cross-dataset run.",
        "Count-based assignment uses only the released counts. LD/AF optimization is the",
        "Saleem attack with its published initialization. SGR is reference haplotype reconstruction.",
        "",
        "## Configuration",
        "",
        f"- HGDP samples: {config['hgdp_sample_count']}",
        f"- beacon sizes: {', '.join(str(value) for value in config['beacon_sizes'])}",
        f"- windows analyzed: {config['windows_analyzed']}",
        f"- replicates: {config['max_replicates']}",
        "",
        "## Carrier F1 and High-Similarity Coverage",
        "",
        "| N | method | n | mean F1 | SD | Sim ≥ 0.95 | mean Hamming | invalid |",
        "|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['beacon_size']} | {row['method_label']} | {row['n_runs']} | {row['mean_f1']:.6f} | "
            f"{row['sd_f1']:.6f} | {row['proportion_similarity_ge_095']:.6f} | {row['mean_hamming']:.2f} | "
            f"{row['invalid_run_count']} |"
        )
    lines.extend(
        [
            "",
            "## Paired Difference, SGR minus Baseline",
            "",
            "| N | baseline | metric | n | mean SGR | mean baseline | mean delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |",
            "|---:|---|---|---:|---:|---:|---:|---|---:|---:|---:|",
        ]
    )
    for row in inference_rows:
        ci = f"[{row['ci95_low']:.6f}, {row['ci95_high']:.6f}]"
        wtl = f"{row['wins']}/{row['ties']}/{row['losses']}"
        lines.append(
            f"| {row['beacon_size']} | {METHOD_LABELS[row['baseline']]} | {row['metric']} | {row['n_pairs']} | "
            f"{row['mean_method']:.6f} | {row['mean_baseline']:.6f} | {row['mean_delta']:.6f} | {ci} | {wtl} | "
            f"{_format_p(float(row['sign_test_p']))} | {_format_p(float(row['wilcoxon_p']))} |"
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "HGDP and 1000 Genomes participants do not overlap. Both cohorts are harmonized",
            "through the gnomAD v3.1.2 HGDP+1KG callset at the retained sites, so this is a",
            "cross-cohort comparison rather than an independent variant-calling comparison.",
            "Windows with fewer than 100 harmonized variants are omitted.",
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def _expected_sgr(path: Path) -> dict[tuple[str, str, int], float]:
    expected: dict[tuple[str, str, int], float] = {}
    if not path.exists():
        return expected
    for row in _read_csv(path):
        if row["method"] != REFERENCE_HAPLOTYPE_RECONSTRUCTION:
            continue
        expected[(row["window_id"], str(row["replicate_id"]), int(row["beacon_size"]))] = float(row["f1"])
    return expected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hgdp-dir", type=Path, default=DEFAULT_HGDP_DIR)
    parser.add_argument("--multiwindow-dir", type=Path, default=DEFAULT_MULTIWINDOW_DIR)
    parser.add_argument("--sgr-metrics", type=Path, default=DEFAULT_SGR_METRICS)
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
    expected_sgr = _expected_sgr(args.sgr_metrics)
    hgdp_sample_count = sum(
        str(row.get("retained", "true")).lower() == "true" for row in _read_csv(args.hgdp_dir / "hgdp_eur_samples.csv")
    )
    rows: list[dict[str, Any]] = []
    windows_analyzed = 0
    for window in windows:
        window_id = window["window_id"]
        if int(overlap_by_window[window_id]["kept_variants"]) < MIN_KEPT_VARIANTS:
            print(f"skip {window_id}", flush=True)
            continue
        windows_analyzed += 1
        hgdp_ids, hgdp_variants, hgdp_positions, hgdp_matrix = _load_hgdp_window(
            args.hgdp_dir / "windows" / f"{window_id}.npz"
        )
        kg_ids, kg_variants, _positions, _alleles, kg_matrix = _load_raw_as_target_carriers(Path(window["raw_path"]))
        kg_index = {variant_id: index for index, variant_id in enumerate(kg_variants)}
        kg_matrix = kg_matrix[:, [kg_index[variant_id] for variant_id in hgdp_variants]]
        for beacon_size in beacon_sizes:
            for replicate_id in replicate_ids:
                seeds = seeds_by_replicate[replicate_id]
                split = split_samples(sample_rows, replicate_id, beacon_size, shared_reference_panel=True)
                reference_ids = split["public_ref"]
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
                        metadata={"window_id": window_id, "experiment": "hgdp_eur_paper_methods"},
                    ),
                )
                truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(beacon_ids))
                runtime_config = {
                    "window_id": window_id,
                    "replicate_id": replicate_id,
                    "beacon_size": beacon_size,
                    "kept_variants": len(hgdp_variants),
                }
                for method in METHOD_ORDER:
                    started = time.perf_counter()
                    result = _with_config(
                        _run_method(method, problem, reference_matrix, beacon_matrix, reference_pairwise, runtime_config),
                        config=runtime_config,
                        optimizer_seed=seeds["optimizer_seed"],
                    )
                    wall_time_sec = time.perf_counter() - started
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
                    if method == REFERENCE_HAPLOTYPE_RECONSTRUCTION:
                        key = (window_id, str(replicate_id), beacon_size)
                        if key in expected_sgr and abs(float(metrics.f1) - expected_sgr[key]) > 1e-9:
                            raise ValueError(
                                f"SGR F1 drifted for {split_id}: {metrics.f1} vs saved {expected_sgr[key]}."
                            )
                    rows.append(
                        {
                            "window_id": window_id,
                            "replicate_id": replicate_id,
                            "beacon_size": beacon_size,
                            "method": method,
                            "method_label": METHOD_LABELS[method],
                            "kept_variants": len(hgdp_variants),
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
                            "wall_time_sec": wall_time_sec,
                        }
                    )
                    print(f"finished {split_id} {method} f1={metrics.f1:.4f} {wall_time_sec:.1f}s", flush=True)

    aggregate_rows = _aggregate(rows)
    inference_rows = _paired_inference(rows)
    config = {
        "hgdp_sample_count": hgdp_sample_count,
        "beacon_sizes": beacon_sizes,
        "windows_analyzed": windows_analyzed,
        "max_replicates": args.max_replicates,
        "methods": list(METHOD_ORDER),
        "saleem_initialization": "official_oracle",
        "saleem_device": "cpu",
    }
    out_dir = args.out_dir
    _write_csv(out_dir / "hgdp_eur_paper_methods_metrics.csv", rows)
    _write_csv(out_dir / "hgdp_eur_paper_methods_aggregate.csv", aggregate_rows)
    _write_csv(out_dir / "hgdp_eur_paper_methods_inference.csv", inference_rows)
    _write_json(out_dir / "hgdp_eur_paper_methods_manifest.json", config)
    _write_report(
        out_dir / "hgdp_eur_paper_methods_report.md",
        aggregate_rows=aggregate_rows,
        inference_rows=inference_rows,
        config=config,
    )


if __name__ == "__main__":
    main()
