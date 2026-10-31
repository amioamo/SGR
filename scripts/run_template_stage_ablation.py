"""Compare reference-only, random-template, and greedy-template ablations.

The full method is Stage 1 greedy template selection followed by Stage 2 exact
count projection. The three baselines remove one or both of those stages:

- reference_only: sample N matched reference rows and do not use y
- random_template_projection: sample those rows, then apply Stage 2
- greedy_template_without_projection: keep the Stage 1 template only
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
from beacon_reconstruction.reconstruction import (
    GREEDY_TEMPLATE_WITHOUT_PROJECTION,
    RANDOM_TEMPLATE_PROJECTION,
    REFERENCE_ONLY,
    REFERENCE_HAPLOTYPE_RECONSTRUCTION,
    TEMPLATE_SELECTION_RULE,
    greedy_template_without_projection,
    random_template_projection,
    reference_only,
    reference_haplotype_reconstruction,
)
from beacon_reconstruction.core import ExperimentConfig, BeaconAttackProblem, BeaconEvaluationTruth, ReconstructionResult
from beacon_reconstruction.evaluation import EvaluationResult, align_rows_hungarian, evaluate_reconstruction
from scripts.experiment_runtime_helpers import attach_reference_panel_mode, run_timed, runtime_metadata_fields
from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers, matrix_for_samples as _matrix_for_samples
from scripts.lib.sample_splits import load_sample_manifest, load_seed_manifest, split_samples
from scripts.lib.reference_helpers import sample_reference_ids as _sample_reference_ids
from scripts.lib.paired_statistics import _summarize_deltas


MANIFEST_ID = "1000g_eur_hg38_template_stage_ablation_v1"
DEFAULT_MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
DEFAULT_OUT_DIR = ROOT / "results" / "block_inheritance_v1" / "template_stage_ablation_shared_panel"
FULL_METHOD = REFERENCE_HAPLOTYPE_RECONSTRUCTION
SIMILARITY_THRESHOLD = 0.95
BOOTSTRAP_SEED = 0
BOOTSTRAP_ITERATIONS = 10_000
METHOD_ORDER = (
    REFERENCE_ONLY,
    RANDOM_TEMPLATE_PROJECTION,
    GREEDY_TEMPLATE_WITHOUT_PROJECTION,
    FULL_METHOD,
)
METHOD_RUNNERS = {
    REFERENCE_ONLY: reference_only,
    RANDOM_TEMPLATE_PROJECTION: random_template_projection,
    GREEDY_TEMPLATE_WITHOUT_PROJECTION: greedy_template_without_projection,
    FULL_METHOD: reference_haplotype_reconstruction,
}


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


def _is_true(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


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


def _record(
    *,
    beacon_size: int,
    window_id: str,
    replicate_id: str,
    method: str,
    panel_size: int,
    result: ReconstructionResult,
    evaluation: EvaluationResult,
    source_overlap_count: int,
    runtime_fields: dict[str, Any],
    similarity: dict[str, float | int],
    projection_edit_fraction_value: float,
    projection_edit_count: int,
) -> dict[str, Any]:
    metrics = evaluation.metrics
    return {
        "beacon_size": beacon_size,
        "window_id": window_id,
        "replicate_id": replicate_id,
        "method": method,
        "method_name": result.method_name,
        "panel_size": panel_size,
        "stage1": result.metadata.get("stage1", TEMPLATE_SELECTION_RULE),
        "stage2": result.metadata.get("stage2", "exact_count_projection"),
        "uses_target_counts": result.metadata.get("uses_target_counts", True),
        "applies_count_projection": result.metadata.get("applies_count_projection", True),
        "source_overlap_count": source_overlap_count,
        "precision": metrics.precision,
        "recall": metrics.recall,
        "f1": metrics.f1,
        "pairwise_error": metrics.pairwise_error,
        "alignment_total_hamming_cost": metrics.alignment_total_hamming_cost,
        "count_mae": metrics.count_mae,
        "max_count_error": metrics.max_count_error,
        "mean_similarity": similarity["mean_similarity"],
        "proportion_similarity_ge_095": similarity["proportion_similarity_ge_095"],
        "n_similarity_ge_095": similarity["n_similarity_ge_095"],
        "n_participants": similarity["n_participants"],
        "projection_edit_fraction": projection_edit_fraction_value,
        "projection_edit_count": projection_edit_count,
        "invalid_run": evaluation.invalid_run.is_invalid,
        "invalid_reasons": ";".join(evaluation.invalid_run.reasons),
        **runtime_fields,
    }


def _mean(rows: list[dict[str, Any]], key: str) -> float:
    return float(np.mean([float(row[key]) for row in rows]))


def projection_edit_count(template: np.ndarray, target_counts: np.ndarray) -> int:
    """Return ``||y - c(T)||_1``, the number of entries Stage 2 must change."""

    return int(np.abs(template.sum(axis=0, dtype=np.int64) - np.asarray(target_counts, dtype=np.int64)).sum())


def projection_edit_fraction(template: np.ndarray, target_counts: np.ndarray) -> float:
    """Return ``||y - c(T)||_1 / (N m)``."""

    row_count, variant_count = template.shape
    if row_count == 0 or variant_count == 0:
        return float("nan")
    return float(projection_edit_count(template, target_counts)) / float(row_count * variant_count)


def participant_similarity_stats(
    true_matrix: np.ndarray,
    reconstructed_matrix: np.ndarray,
    *,
    threshold: float = SIMILARITY_THRESHOLD,
) -> dict[str, float | int]:
    """Participant similarity after the same Hungarian alignment used for carrier F1."""

    alignment = align_rows_hungarian(true_matrix, reconstructed_matrix)
    variant_count = int(true_matrix.shape[1])
    if variant_count == 0:
        return {
            "mean_similarity": float("nan"),
            "proportion_similarity_ge_095": float("nan"),
            "n_similarity_ge_095": 0,
            "n_participants": int(true_matrix.shape[0]),
        }
    true_aligned = true_matrix[alignment.true_indices]
    reconstructed_aligned = reconstructed_matrix[alignment.reconstructed_indices]
    hamming = np.not_equal(true_aligned, reconstructed_aligned).sum(axis=1)
    similarity = 1.0 - hamming.astype(float) / float(variant_count)
    at_threshold = similarity >= threshold
    return {
        "mean_similarity": float(np.mean(similarity)),
        "proportion_similarity_ge_095": float(np.mean(at_threshold)),
        "n_similarity_ge_095": int(np.sum(at_threshold)),
        "n_participants": int(true_matrix.shape[0]),
    }


def _aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    keys = sorted(
        {(int(row["beacon_size"]), row["method"]) for row in rows},
        key=lambda item: (item[0], METHOD_ORDER.index(item[1])),
    )
    for beacon_size, method in keys:
        selected = [row for row in rows if int(row["beacon_size"]) == beacon_size and row["method"] == method]
        out.append(
            {
                "beacon_size": beacon_size,
                "method": method,
                "panel_size": selected[0]["panel_size"],
                "n_runs": len(selected),
                "mean_f1": _mean(selected, "f1"),
                "median_f1": float(np.median([float(row["f1"]) for row in selected])),
                "mean_hamming": _mean(selected, "alignment_total_hamming_cost"),
                "mean_count_mae": _mean(selected, "count_mae"),
                "mean_max_count_error": _mean(selected, "max_count_error"),
                "mean_similarity": _mean(selected, "mean_similarity"),
                "proportion_similarity_ge_095": _mean(selected, "proportion_similarity_ge_095"),
                "mean_projection_edit_fraction": _mean(selected, "projection_edit_fraction"),
                "mean_wall_time_sec": _mean(selected, "wall_time_sec"),
                "total_source_overlap_count": int(sum(int(row["source_overlap_count"]) for row in selected)),
                "invalid_run_count": int(sum(_is_true(row["invalid_run"]) for row in selected)),
            }
        )
    return out


def _deltas(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {
        (int(row["beacon_size"]), row["window_id"], row["replicate_id"], row["method"]): row for row in rows
    }
    units = sorted({(int(row["beacon_size"]), row["window_id"], row["replicate_id"]) for row in rows})
    out: list[dict[str, Any]] = []
    for beacon_size, window_id, replicate_id in units:
        full = by_key[(beacon_size, window_id, replicate_id, FULL_METHOD)]
        for method in METHOD_ORDER:
            if method == FULL_METHOD:
                continue
            row = by_key[(beacon_size, window_id, replicate_id, method)]
            out.append(
                {
                    "beacon_size": beacon_size,
                    "window_id": window_id,
                    "replicate_id": replicate_id,
                    "method": method,
                    "delta_f1_vs_full": float(row["f1"]) - float(full["f1"]),
                    "delta_hamming_vs_full": float(row["alignment_total_hamming_cost"])
                    - float(full["alignment_total_hamming_cost"]),
                }
            )
    return out


def _paired_inference(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Paired tests within each beacon size. Delta is primary minus comparator."""

    by_key = {
        (int(row["beacon_size"]), row["window_id"], row["replicate_id"], row["method"]): row for row in rows
    }
    beacon_sizes = sorted({int(row["beacon_size"]) for row in rows})
    summaries: list[dict[str, Any]] = []
    carrier_metrics = (
        ("f1", True),
        ("proportion_similarity_ge_095", True),
    )
    for beacon_size in beacon_sizes:
        units = sorted(
            {
                (row["window_id"], row["replicate_id"])
                for row in rows
                if int(row["beacon_size"]) == beacon_size
            }
        )
        for metric, higher_is_better in carrier_metrics:
            for baseline in METHOD_ORDER:
                if baseline == FULL_METHOD:
                    continue
                deltas: list[float] = []
                baseline_values: list[float] = []
                method_values: list[float] = []
                for window_id, replicate_id in units:
                    full = by_key[(beacon_size, window_id, replicate_id, FULL_METHOD)]
                    other = by_key[(beacon_size, window_id, replicate_id, baseline)]
                    method_values.append(float(full[metric]))
                    baseline_values.append(float(other[metric]))
                    deltas.append(float(full[metric]) - float(other[metric]))
                summary = _summarize_deltas(
                    comparison=f"{FULL_METHOD}_minus_{baseline}",
                    metric=metric,
                    deltas=deltas,
                    baseline_values=baseline_values,
                    method_values=method_values,
                    higher_is_better=higher_is_better,
                    bootstrap_seed=BOOTSTRAP_SEED,
                    bootstrap_iterations=BOOTSTRAP_ITERATIONS,
                )
                summary["beacon_size"] = beacon_size
                summaries.append(summary)

        edit_deltas: list[float] = []
        random_values: list[float] = []
        greedy_values: list[float] = []
        for window_id, replicate_id in units:
            random_row = by_key[(beacon_size, window_id, replicate_id, REFERENCE_ONLY)]
            greedy_row = by_key[(beacon_size, window_id, replicate_id, GREEDY_TEMPLATE_WITHOUT_PROJECTION)]
            random_values.append(float(random_row["projection_edit_fraction"]))
            greedy_values.append(float(greedy_row["projection_edit_fraction"]))
            edit_deltas.append(float(greedy_row["projection_edit_fraction"]) - float(random_row["projection_edit_fraction"]))
        edit_summary = _summarize_deltas(
            comparison=f"{GREEDY_TEMPLATE_WITHOUT_PROJECTION}_minus_{REFERENCE_ONLY}",
            metric="projection_edit_fraction",
            deltas=edit_deltas,
            baseline_values=random_values,
            method_values=greedy_values,
            higher_is_better=False,
            bootstrap_seed=BOOTSTRAP_SEED,
            bootstrap_iterations=BOOTSTRAP_ITERATIONS,
        )
        edit_summary["beacon_size"] = beacon_size
        summaries.append(edit_summary)
    return summaries


def _format_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.4g}"


def _write_report(
    path: Path,
    *,
    aggregate_rows: list[dict[str, Any]],
    delta_rows: list[dict[str, Any]],
    inference_rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> None:
    lines = [
        "# Stage Ablation On The Matched Reference Panel",
        "",
        "## Configuration",
        "",
        f"- windows: {config['max_windows']}",
        f"- replicates: {config['max_replicates']}",
        f"- beacon sizes: {','.join(str(size) for size in config['beacon_sizes'])}",
        f"- panel size: {config['panel_size']}",
        f"- reference panel mode: {config['reference_panel_mode']}",
        "",
        "Stage 1 of the full method is greedy template selection. Stage 2 is exact count projection.",
        "`reference_only` samples N panel rows and does not use y.",
        "`random_template_projection` uses that same random draw as the template, then applies Stage 2.",
        "`greedy_template_without_projection` keeps the Stage 1 template and skips Stage 2.",
        "Count violations are expected for the two methods that skip projection; their F1 is still reported.",
        "",
        "## Aggregate Metrics",
        "",
        "| N | method | n_runs | mean_f1 | median_f1 | mean_hamming | mean_count_mae | mean_max_count_error | invalid |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['beacon_size']} | {row['method']} | {row['n_runs']} | "
            f"{row['mean_f1']:.6f} | {row['median_f1']:.6f} | {row['mean_hamming']:.2f} | "
            f"{row['mean_count_mae']:.4f} | {row['mean_max_count_error']:.4f} | {row['invalid_run_count']} |"
        )

    lines.extend(
        [
            "",
            "## Gain Versus Each Baseline",
            "",
            "Positive `full_minus_baseline_f1` is the F1 gained by using both stages instead of that baseline.",
            "`greedy_template_without_projection` isolates the Stage 2 projection gain.",
            "`random_template_projection` isolates the Stage 1 greedy-selection gain.",
            "",
            "| N | baseline | mean_baseline_f1 | mean_full_f1 | full_minus_baseline_f1 | wins_for_full | ties | losses |",
            "|---:|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for beacon_size, method in sorted(
        {(int(row["beacon_size"]), row["method"]) for row in delta_rows},
        key=lambda item: (item[0], METHOD_ORDER.index(item[1])),
    ):
        selected = [row for row in delta_rows if int(row["beacon_size"]) == beacon_size and row["method"] == method]
        deltas = np.asarray([float(row["delta_f1_vs_full"]) for row in selected], dtype=float)
        full_f1 = next(
            float(row["mean_f1"])
            for row in aggregate_rows
            if int(row["beacon_size"]) == beacon_size and row["method"] == FULL_METHOD
        )
        baseline_f1 = next(
            float(row["mean_f1"])
            for row in aggregate_rows
            if int(row["beacon_size"]) == beacon_size and row["method"] == method
        )
        gain = full_f1 - baseline_f1
        lines.append(
            f"| {beacon_size} | {method} | {baseline_f1:.6f} | {full_f1:.6f} | {gain:.6f} | "
            f"{int(np.sum(deltas < 0))} | {int(np.sum(deltas == 0))} | {int(np.sum(deltas > 0))} |"
        )

    lines.extend(
        [
            "",
            "## Participant Similarity",
            "",
            "Similarity is `1 - Hamming / m` after the same one-to-one Hungarian alignment used for carrier F1.",
            f"`sim>=0.95` is the fraction of aligned participants at or above {SIMILARITY_THRESHOLD:.2f}.",
            "",
            "| N | method | n_runs | participants | mean_similarity | sim>=0.95 |",
            "|---:|---|---:|---:|---:|---:|",
        ]
    )
    for row in aggregate_rows:
        participants = int(row["n_runs"]) * int(row["beacon_size"])
        lines.append(
            f"| {row['beacon_size']} | {row['method']} | {row['n_runs']} | {participants} | "
            f"{row['mean_similarity']:.6f} | {row['proportion_similarity_ge_095']:.6f} |"
        )

    lines.extend(
        [
            "",
            "## Paired Inference",
            "",
            "Each pair is one window by replicate. "
            f"The 95% CI is a percentile bootstrap of the mean paired difference "
            f"({BOOTSTRAP_ITERATIONS} resamples, seed {BOOTSTRAP_SEED}). "
            "The sign test is two-sided and drops ties. "
            "The Wilcoxon signed-rank test is two-sided and drops exact zeros.",
            "For F1 and sim>=0.95, delta is full method minus baseline, and a win means the full method is higher.",
            "Count-violation runs are retained.",
            "",
            "| N | comparison | metric | n | mean_delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |",
            "|---:|---|---|---:|---:|---|---:|---:|---:|",
        ]
    )
    inference_order = {"f1": 0, "proportion_similarity_ge_095": 1, "projection_edit_fraction": 2}
    for row in sorted(
        inference_rows,
        key=lambda item: (
            int(item["beacon_size"]),
            inference_order.get(item["metric"], 9),
            item["comparison"],
        ),
    ):
        if row["metric"] == "projection_edit_fraction":
            continue
        ci = f"[{row['ci95_low']:.6f}, {row['ci95_high']:.6f}]"
        wtl = f"{row['wins']}/{row['ties']}/{row['losses']}"
        lines.append(
            f"| {row['beacon_size']} | {row['comparison']} | {row['metric']} | {row['n_pairs']} | "
            f"{row['mean_delta']:.6f} | {ci} | {wtl} | {_format_p(float(row['sign_test_p']))} | "
            f"{_format_p(float(row['wilcoxon_p']))} |"
        )

    lines.extend(
        [
            "",
            "## Projection Edit Fraction",
            "",
            "For a template `T`, the edit fraction is `||y - c(T)||_1 / (N m)`.",
            "Stage 2 changes exactly that fraction of entries, because each column flips `|y_j - c_j(T)|` bits.",
            "The random template is the `reference_only` sample. The greedy template is Stage 1.",
            "Delta is greedy minus random, so a negative value means the greedy template needs fewer edits.",
            "A win means the greedy template has the smaller edit fraction.",
            "",
            "| N | template | mean_edit_fraction |",
            "|---:|---|---:|",
        ]
    )
    for row in aggregate_rows:
        if row["method"] not in {REFERENCE_ONLY, GREEDY_TEMPLATE_WITHOUT_PROJECTION}:
            continue
        lines.append(
            f"| {row['beacon_size']} | {row['method']} | {row['mean_projection_edit_fraction']:.6f} |"
        )
    lines.extend(
        [
            "",
            "| N | comparison | n | mean_delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |",
            "|---:|---|---:|---:|---|---:|---:|---:|",
        ]
    )
    for row in inference_rows:
        if row["metric"] != "projection_edit_fraction":
            continue
        ci = f"[{row['ci95_low']:.6f}, {row['ci95_high']:.6f}]"
        wtl = f"{row['wins']}/{row['ties']}/{row['losses']}"
        lines.append(
            f"| {row['beacon_size']} | {row['comparison']} | {row['n_pairs']} | {row['mean_delta']:.6f} | "
            f"{ci} | {wtl} | {_format_p(float(row['sign_test_p']))} | {_format_p(float(row['wilcoxon_p']))} |"
        )

    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Rows are aligned with Hungarian matching before F1 is computed, including runs whose column sums",
            "do not match y. `invalid_run` records that count violation and does not drop the F1.",
            f"Reference panel mode: `{config['reference_panel_mode']}`.",
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--beacon-sizes", default="25,50,100,250")
    parser.add_argument("--panel-size", type=int, default=100)
    parser.add_argument("--shared-reference-panel", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--multiwindow-dir", type=Path, default=DEFAULT_MULTIWINDOW_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    out_dir = args.out_dir
    multiwindow_dir = args.multiwindow_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    windows = _read_csv(multiwindow_dir / "multiwindow_manifest.csv")[: max(0, args.max_windows)]
    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]
    beacon_sizes = [int(value) for value in args.beacon_sizes.split(",") if value.strip()]
    config = {
        "max_windows": args.max_windows,
        "max_replicates": args.max_replicates,
        "beacon_sizes": beacon_sizes,
        "panel_size": args.panel_size,
        "reference_panel_mode": "shared" if args.shared_reference_panel else "dual_disjoint",
        "template_selection_rule": TEMPLATE_SELECTION_RULE,
        "methods": list(METHOD_ORDER),
        "similarity_threshold": SIMILARITY_THRESHOLD,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
        "multiwindow_dir": str(multiwindow_dir),
        "out_dir": str(out_dir),
    }

    rows: list[dict[str, Any]] = []
    for window in windows:
        window_id = window["window_id"]
        sample_ids, variant_ids, positions, _target_alleles, carrier_matrix = _load_raw_as_target_carriers(
            Path(window["raw_path"])
        )
        for beacon_size in beacon_sizes:
            for replicate_id in replicate_ids:
                seeds = seeds_by_replicate[replicate_id]
                split = split_samples(
                    sample_rows,
                    replicate_id,
                    beacon_size,
                    shared_reference_panel=args.shared_reference_panel,
                )
                template_ref_ids = split["public_ref"] if args.shared_reference_panel else split["arg_ref"]
                panel_ids = _sample_reference_ids(
                    template_ref_ids,
                    panel_size=args.panel_size,
                    seed=seeds["group_selection_seed"] + args.panel_size,
                )
                panel_matrix = _matrix_for_samples(sample_ids, carrier_matrix, panel_ids)
                beacon_matrix = _matrix_for_samples(sample_ids, carrier_matrix, split["beacon"])
                public_ref_matrix = _matrix_for_samples(sample_ids, carrier_matrix, split["public_ref"])
                target_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                reference_pairwise = pairwise_moments(public_ref_matrix)
                split_id = f"{window_id}_replicate_{replicate_id}_N{beacon_size}"
                runtime_config = {
                    "window_id": window_id,
                    "replicate_id": replicate_id,
                    "beacon_size": beacon_size,
                    "panel_size": args.panel_size,
                    "config": config,
                }
                problem = BeaconAttackProblem(
                    manifest_id=MANIFEST_ID,
                    split_id=split_id,
                    beacon_size=beacon_size,
                    variant_ids=variant_ids,
                    positions=positions,
                    query_responses=target_counts,
                    target_counts=target_counts,
                    public_ld=reference_pairwise,
                    config=ExperimentConfig(
                        initialization_seed=seeds["initialization_seed"],
                        optimizer_seed=seeds["optimizer_seed"],
                        group_selection_seed=seeds["group_selection_seed"],
                        metadata={"window_id": window_id, "experiment": "template_stage_ablation"},
                    ),
                )
                truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(split["beacon"]))
                overlap_count = len(set(split["beacon"]).intersection(panel_ids))
                completed: dict[str, tuple[ReconstructionResult, EvaluationResult, float]] = {}
                for method in METHOD_ORDER:
                    runner = METHOD_RUNNERS[method]
                    result, wall_time_sec = run_timed(
                        lambda runner=runner, problem=problem, panel_matrix=panel_matrix, runtime_config=runtime_config: runner(
                            problem,
                            panel_matrix,
                            metadata=runtime_config,
                        )
                    )
                    result = _with_config(result, config=runtime_config, optimizer_seed=seeds["optimizer_seed"])
                    result = attach_reference_panel_mode(result, shared=args.shared_reference_panel)
                    evaluation = evaluate_reconstruction(
                        result,
                        truth,
                        target_counts=target_counts,
                        reference_pairwise=reference_pairwise,
                        source_sets={"beacon": split["beacon"], "template_reference": panel_ids},
                        accepted_manifest_ids=None,
                    )
                    completed[method] = (result, evaluation, wall_time_sec)

                random_template = completed[REFERENCE_ONLY][0].reconstructed_binary
                greedy_template = completed[GREEDY_TEMPLATE_WITHOUT_PROJECTION][0].reconstructed_binary
                random_projected = completed[RANDOM_TEMPLATE_PROJECTION][0].reconstructed_binary
                greedy_projected = completed[FULL_METHOD][0].reconstructed_binary
                random_edit_count = projection_edit_count(random_template, target_counts)
                greedy_edit_count = projection_edit_count(greedy_template, target_counts)
                random_edit = random_edit_count / float(beacon_matrix.size)
                greedy_edit = greedy_edit_count / float(beacon_matrix.size)
                if random_edit_count != int(np.not_equal(random_template, random_projected).sum()):
                    raise ValueError(f"Random-template edit count mismatch for {split_id}.")
                if greedy_edit_count != int(np.not_equal(greedy_template, greedy_projected).sum()):
                    raise ValueError(f"Greedy-template edit count mismatch for {split_id}.")
                edit_fraction_by_method = {
                    REFERENCE_ONLY: random_edit,
                    RANDOM_TEMPLATE_PROJECTION: random_edit,
                    GREEDY_TEMPLATE_WITHOUT_PROJECTION: greedy_edit,
                    FULL_METHOD: greedy_edit,
                }
                edit_count_by_method = {
                    REFERENCE_ONLY: random_edit_count,
                    RANDOM_TEMPLATE_PROJECTION: random_edit_count,
                    GREEDY_TEMPLATE_WITHOUT_PROJECTION: greedy_edit_count,
                    FULL_METHOD: greedy_edit_count,
                }
                for method, (result, evaluation, wall_time_sec) in completed.items():
                    rows.append(
                        _record(
                            beacon_size=beacon_size,
                            window_id=window_id,
                            replicate_id=replicate_id,
                            method=method,
                            panel_size=len(panel_ids),
                            result=result,
                            evaluation=evaluation,
                            source_overlap_count=overlap_count,
                            runtime_fields=runtime_metadata_fields(result, wall_time_sec),
                            similarity=participant_similarity_stats(beacon_matrix, result.reconstructed_binary),
                            projection_edit_fraction_value=edit_fraction_by_method[method],
                            projection_edit_count=edit_count_by_method[method],
                        )
                    )

    aggregate_rows = _aggregate(rows)
    delta_rows = _deltas(rows)
    inference_rows = _paired_inference(rows)
    _write_csv(out_dir / "template_stage_ablation_metrics.csv", rows)
    _write_csv(out_dir / "template_stage_ablation_aggregate.csv", aggregate_rows)
    _write_csv(out_dir / "template_stage_ablation_deltas.csv", delta_rows)
    _write_csv(out_dir / "template_stage_ablation_inference.csv", inference_rows)
    _write_json(out_dir / "template_stage_ablation_manifest.json", config)
    _write_report(
        out_dir / "template_stage_ablation_report.md",
        aggregate_rows=aggregate_rows,
        delta_rows=delta_rows,
        inference_rows=inference_rows,
        config=config,
    )


if __name__ == "__main__":
    main()
