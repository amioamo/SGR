"""Nested reference-panel-size sweep for reference haplotype reconstruction.

Within each window, replicate, and beacon size, the beacon is fixed. The
reference is a prefix of one permutation of the beacon's complement, so a
larger panel contains every smaller panel. This is separate from
``reference_panel_size_gradient``, which reserves an extra public panel of 100
before taking the template prefix and also runs pairwise refinement.
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
from beacon_reconstruction.reconstruction import reference_haplotype_reconstruction
from beacon_reconstruction.core import ExperimentConfig, BeaconAttackProblem, BeaconEvaluationTruth, ReconstructionResult
from beacon_reconstruction.evaluation import EvaluationResult, evaluate_reconstruction
from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers
from scripts.run_multiancestry_reference_matching import _matrix_for_indices
from scripts.run_template_stage_ablation import participant_similarity_stats
from scripts.lib.paired_statistics import _summarize_deltas


MANIFEST_ID = "1000g_eur_hg38_nested_reference_panel_size_v1"
DEFAULT_PANEL_DIR = ROOT / "data" / "multiwindow" / "ancestry_panels"
DEFAULT_OUT_DIR = ROOT / "results" / "block_inheritance_v1" / "reference_panel_size_nested_eur"
BASELINE_PANEL_SIZE = 100
BOOTSTRAP_SEED = 0
BOOTSTRAP_ITERATIONS = 10_000


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _parse_int_list(raw: str) -> list[int]:
    return sorted({int(value.strip()) for value in raw.split(",") if value.strip()})


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


def _nested_split(total: int, *, beacon_size: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    if beacon_size <= 0 or beacon_size >= total:
        raise ValueError(f"beacon_size={beacon_size} is not a valid split of {total} samples.")
    rng = np.random.default_rng(seed)
    selected = rng.permutation(total).astype(int)
    beacon = np.sort(selected[:beacon_size])
    complement = selected[beacon_size:]
    return beacon, complement


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
    keys = sorted({(int(row["beacon_size"]), int(row["panel_size"])) for row in rows})
    out: list[dict[str, Any]] = []
    for beacon_size, panel_size in keys:
        selected = [
            row
            for row in rows
            if int(row["beacon_size"]) == beacon_size and int(row["panel_size"]) == panel_size
        ]
        out.append(
            {
                "beacon_size": beacon_size,
                "panel_size": panel_size,
                "n_runs": len(selected),
                "mean_f1": _mean(selected, "f1"),
                "median_f1": float(np.median([float(row["f1"]) for row in selected])),
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
        (int(row["beacon_size"]), int(row["panel_size"]), row["window_id"], row["replicate_id"]): row for row in rows
    }
    panel_sizes = sorted({int(row["panel_size"]) for row in rows})
    beacon_sizes = sorted({int(row["beacon_size"]) for row in rows})
    units = sorted({(row["window_id"], row["replicate_id"]) for row in rows})
    inference: list[dict[str, Any]] = []
    for beacon_size in beacon_sizes:
        for panel_size in panel_sizes:
            if panel_size == BASELINE_PANEL_SIZE:
                continue
            f1_deltas: list[float] = []
            f1_base: list[float] = []
            f1_method: list[float] = []
            sim_deltas: list[float] = []
            sim_base: list[float] = []
            sim_method: list[float] = []
            for window_id, replicate_id in units:
                baseline = by_key.get((beacon_size, BASELINE_PANEL_SIZE, window_id, replicate_id))
                method = by_key.get((beacon_size, panel_size, window_id, replicate_id))
                if baseline is None or method is None:
                    continue
                f1_deltas.append(float(method["f1"]) - float(baseline["f1"]))
                f1_base.append(float(baseline["f1"]))
                f1_method.append(float(method["f1"]))
                sim_deltas.append(
                    float(method["proportion_similarity_ge_095"]) - float(baseline["proportion_similarity_ge_095"])
                )
                sim_base.append(float(baseline["proportion_similarity_ge_095"]))
                sim_method.append(float(method["proportion_similarity_ge_095"]))
            comparison = f"R{panel_size}_minus_R{BASELINE_PANEL_SIZE}"
            seed = BOOTSTRAP_SEED + beacon_size * 1000 + panel_size
            for metric, deltas, base, method_values in (
                ("f1", f1_deltas, f1_base, f1_method),
                ("proportion_similarity_ge_095", sim_deltas, sim_base, sim_method),
            ):
                summary = _summarize_deltas(
                    comparison=comparison,
                    metric=metric,
                    deltas=deltas,
                    baseline_values=base,
                    method_values=method_values,
                    higher_is_better=True,
                    bootstrap_seed=seed,
                    bootstrap_iterations=BOOTSTRAP_ITERATIONS,
                )
                summary["beacon_size"] = beacon_size
                summary["panel_size"] = panel_size
                summary["baseline_panel_size"] = BASELINE_PANEL_SIZE
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
        "# Nested Reference Panel Size",
        "",
        "The beacon is drawn once per window, replicate, and beacon size. The reference",
        "is a nested prefix of the complement, so every larger panel contains the smaller",
        "panels. The method is reference haplotype reconstruction (SGR). Deltas are F1(R) minus F1(R=100)",
        "on the same beacon.",
        "",
        "## Configuration",
        "",
        f"- ancestry: {config['ancestry']}",
        f"- beacon sizes: {', '.join(str(value) for value in config['beacon_sizes'])}",
        f"- panel sizes: {', '.join(str(value) for value in config['panel_sizes'])}",
        f"- windows: {config['max_windows']}",
        f"- replicates: {config['replicates']}",
        "",
        "## Carrier F1",
        "",
        "| N | R | n | mean F1 | median F1 | mean Hamming | Sim ≥ 0.95 | overlap | invalid |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['beacon_size']} | {row['panel_size']} | {row['n_runs']} | {row['mean_f1']:.6f} | "
            f"{row['median_f1']:.6f} | {row['mean_hamming']:.2f} | {row['proportion_similarity_ge_095']:.6f} | "
            f"{row['total_source_overlap_count']} | {row['invalid_run_count']} |"
        )
    lines.extend(
        [
            "",
            "## Paired Difference versus R = 100",
            "",
            "Positive means the listed panel size has the higher value. Intervals are percentile",
            "bootstrap 95% CIs of the mean paired difference.",
            "",
            "| N | R | metric | n | mean delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |",
            "|---:|---:|---|---:|---:|---|---:|---:|---:|",
        ]
    )
    for row in inference_rows:
        ci = f"[{row['ci95_low']:.6f}, {row['ci95_high']:.6f}]"
        wtl = f"{row['wins']}/{row['ties']}/{row['losses']}"
        lines.append(
            f"| {row['beacon_size']} | {row['panel_size']} | {row['metric']} | {row['n_pairs']} | "
            f"{row['mean_delta']:.6f} | {ci} | {wtl} | {_format_p(float(row['sign_test_p']))} | "
            f"{_format_p(float(row['wilcoxon_p']))} |"
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "The beacon and the reference are disjoint samples from the same 1000 Genomes EUR panel.",
            "This sweep varies reference size only. It is not a cross-dataset result.",
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-dir", type=Path, default=DEFAULT_PANEL_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--ancestry", default="EUR")
    parser.add_argument("--beacon-sizes", default="25,50,100")
    parser.add_argument("--panel-sizes", default="10,25,50,100,200,400")
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--replicates", type=int, default=3)
    args = parser.parse_args()

    ancestry = args.ancestry.upper()
    beacon_sizes = _parse_int_list(args.beacon_sizes)
    panel_sizes = _parse_int_list(args.panel_sizes)
    if BASELINE_PANEL_SIZE not in panel_sizes:
        raise ValueError(f"Panel sizes must include the baseline R={BASELINE_PANEL_SIZE}.")
    panel_manifest = _read_csv(args.panel_dir / "ancestry_panel_manifest.csv")
    windows = [row for row in panel_manifest if row["ancestry"].upper() == ancestry][: max(0, args.max_windows)]
    if len(windows) != args.max_windows:
        raise ValueError(f"Expected {args.max_windows} {ancestry} windows, found {len(windows)}.")

    config = {
        "ancestry": ancestry,
        "beacon_sizes": beacon_sizes,
        "panel_sizes": panel_sizes,
        "baseline_panel_size": BASELINE_PANEL_SIZE,
        "max_windows": args.max_windows,
        "replicates": args.replicates,
        "panel_dir": str(args.panel_dir),
        "out_dir": str(args.out_dir),
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
    }
    rows: list[dict[str, Any]] = []
    for window in windows:
        window_id = window["window_id"]
        sample_ids, variant_ids, positions, _target_alleles, carrier_matrix = _load_raw_as_target_carriers(
            Path(window["raw_path"])
        )
        for beacon_size in beacon_sizes:
            for replicate_id in range(args.replicates):
                seed_base = 100_000 + 10_000 * replicate_id + 100 * beacon_size + len(window_id)
                beacon_indices, complement = _nested_split(len(sample_ids), beacon_size=beacon_size, seed=seed_base)
                if max(panel_sizes) > len(complement):
                    raise ValueError(
                        f"R={max(panel_sizes)} exceeds the complement of N={beacon_size} "
                        f"({len(complement)} samples) in {window_id}."
                    )
                beacon_matrix = _matrix_for_indices(carrier_matrix, beacon_indices)
                target_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                beacon_sample_ids = [sample_ids[int(index)] for index in beacon_indices]
                split_id = f"{window_id}_{ancestry}_replicate_{replicate_id}_N{beacon_size}"
                problem = BeaconAttackProblem(
                    manifest_id=MANIFEST_ID,
                    split_id=split_id,
                    beacon_size=beacon_size,
                    variant_ids=variant_ids,
                    positions=positions,
                    query_responses=target_counts,
                    target_counts=target_counts,
                    public_ld=pairwise_moments(_matrix_for_indices(carrier_matrix, complement[:BASELINE_PANEL_SIZE])),
                    config=ExperimentConfig(
                        initialization_seed=seed_base + 1,
                        optimizer_seed=seed_base + 2,
                        group_selection_seed=seed_base + 3,
                        metadata={"window_id": window_id, "experiment": "nested_reference_panel_size"},
                    ),
                )
                truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(beacon_sample_ids))
                for panel_size in panel_sizes:
                    reference_indices = np.sort(complement[:panel_size].astype(int))
                    reference_ids = [sample_ids[int(index)] for index in reference_indices]
                    overlap = len(set(beacon_sample_ids).intersection(reference_ids))
                    if overlap:
                        raise ValueError(f"Beacon/reference overlap at {split_id}, R={panel_size}.")
                    reference_matrix = _matrix_for_indices(carrier_matrix, reference_indices)
                    runtime_config = {
                        "window_id": window_id,
                        "beacon_size": beacon_size,
                        "panel_size": panel_size,
                        "replicate_id": replicate_id,
                    }
                    result = _with_config(
                        reference_haplotype_reconstruction(problem, reference_matrix, metadata=runtime_config),
                        config=runtime_config,
                        optimizer_seed=seed_base + 2,
                    )
                    evaluation: EvaluationResult = evaluate_reconstruction(
                        result,
                        truth,
                        target_counts=target_counts,
                        reference_pairwise=problem.public_ld,
                        source_sets={"beacon": beacon_sample_ids, "template_reference": reference_ids},
                        accepted_manifest_ids=None,
                    )
                    similarity = participant_similarity_stats(beacon_matrix, result.reconstructed_binary)
                    metrics = evaluation.metrics
                    rows.append(
                        {
                            "window_id": window_id,
                            "replicate_id": replicate_id,
                            "beacon_size": beacon_size,
                            "panel_size": panel_size,
                            "method": result.method_name,
                            "complement_size": int(len(complement)),
                            "precision": metrics.precision,
                            "recall": metrics.recall,
                            "f1": metrics.f1,
                            "alignment_total_hamming_cost": metrics.alignment_total_hamming_cost,
                            "count_mae": metrics.count_mae,
                            "max_count_error": metrics.max_count_error,
                            "mean_similarity": similarity["mean_similarity"],
                            "proportion_similarity_ge_095": similarity["proportion_similarity_ge_095"],
                            "source_overlap_count": overlap,
                            "invalid_run": evaluation.invalid_run.is_invalid,
                            "invalid_reasons": ";".join(evaluation.invalid_run.reasons),
                        }
                    )
                print(f"finished {split_id}", flush=True)

    aggregate_rows = _aggregate(rows)
    inference_rows = _paired_inference(rows)
    out_dir = args.out_dir
    _write_csv(out_dir / "nested_reference_panel_size_metrics.csv", rows)
    _write_csv(out_dir / "nested_reference_panel_size_aggregate.csv", aggregate_rows)
    _write_csv(out_dir / "nested_reference_panel_size_inference.csv", inference_rows)
    _write_json(out_dir / "nested_reference_panel_size_manifest.json", config)
    _write_report(
        out_dir / "nested_reference_panel_size_report.md",
        aggregate_rows=aggregate_rows,
        inference_rows=inference_rows,
        config=config,
    )


if __name__ == "__main__":
    main()
