"""Compute per-participant privacy metrics for genomic data-sharing reconstruction."""

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

from beacon_reconstruction.baselines import independent_count, pairwise_moments, pairwise_refinement
from beacon_reconstruction.reconstruction import reference_haplotype_reconstruction
from beacon_reconstruction.core import ExperimentConfig, BeaconAttackProblem, BeaconEvaluationTruth, ReconstructionResult
from beacon_reconstruction.evaluation import align_rows_hungarian, compute_count_metrics, evaluate_reconstruction
from beacon_reconstruction.ld_af_baseline import SALEEM_RECONSTRUCTION_ATTACK
from scripts.ld_af_helpers import (
    add_saleem_arguments,
    merge_saleem_rows,
    run_saleem_attack,
    saleem_config_dict,
)
from scripts.lib.multiwindow_data import (
    load_raw_carriers,
    matrix_for_samples,
)
from scripts.lib.sample_splits import load_sample_manifest, load_seed_manifest, split_samples
from scripts.lib.reference_helpers import (
    sample_reference_ids,
    shuffle_reference_columns,
)


MANIFEST_ID = "1000g_eur_hg38_participant_privacy_v1"
MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
OUT_DIR = ROOT / "results" / "experiments" / "row_privacy_template_count"
DEFAULT_SHARED_OUT_DIR = ROOT / "results" / "experiments" / "row_privacy_template_count_shared_panel"
SIMILARITY_THRESHOLDS = (0.90, 0.95, 0.99, 1.0)


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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_serialize_value(payload), indent=2, sort_keys=True) + "\n")


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


def _safe_divide(numerator: int | float, denominator: int | float) -> float:
    if denominator == 0:
        return float("nan")
    return float(numerator / denominator)


def _row_f1(true_row: np.ndarray, reconstructed_row: np.ndarray) -> tuple[float, float, float, int, int, int]:
    true_binary = np.asarray(true_row, dtype=np.int8)
    reconstructed_binary = np.asarray(reconstructed_row, dtype=np.int8)
    tp = int(np.logical_and(true_binary == 1, reconstructed_binary == 1).sum())
    fp = int(np.logical_and(true_binary == 0, reconstructed_binary == 1).sum())
    fn = int(np.logical_and(true_binary == 1, reconstructed_binary == 0).sum())
    precision = _safe_divide(tp, tp + fp)
    recall = _safe_divide(tp, tp + fn)
    if np.isfinite(precision) and np.isfinite(recall) and (precision + recall) > 0:
        f1 = float(2 * precision * recall / (precision + recall))
    else:
        f1 = float("nan")
    return precision, recall, f1, tp, fp, fn


def _row_metric_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        int(row["beacon_size"]),
        row["window_id"],
        row["replicate_id"],
        row["method"],
        row["beacon_sample_id"],
    )


def _run_summary_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        int(row["beacon_size"]),
        row["window_id"],
        row["replicate_id"],
        row["method"],
    )


def _merge_existing_rows(
    new_rows: list[dict[str, Any]],
    *,
    metrics_path: Path,
    key_fn: Any,
) -> list[dict[str, Any]]:
    if not metrics_path.exists():
        return new_rows
    existing_rows = _read_csv(metrics_path)
    new_keys = {key_fn(row) for row in new_rows}
    kept = [row for row in existing_rows if key_fn(row) not in new_keys]
    return kept + new_rows


def _row_privacy_rows(
    *,
    beacon_size: int,
    window_id: str,
    replicate_id: str,
    method: str,
    reference_population: str,
    template_design: str,
    panel_size: int | None,
    beacon_ids: list[str],
    true_matrix: np.ndarray,
    reconstructed_matrix: np.ndarray,
    max_count_error: float,
    source_overlap_count: int,
    invalid_run: bool,
) -> list[dict[str, Any]]:
    alignment = align_rows_hungarian(true_matrix, reconstructed_matrix)
    out: list[dict[str, Any]] = []
    variant_count = int(true_matrix.shape[1])
    for true_index, reconstructed_index in zip(alignment.true_indices, alignment.reconstructed_indices):
        true_row = true_matrix[int(true_index)]
        reconstructed_row = reconstructed_matrix[int(reconstructed_index)]
        hamming = int(np.not_equal(true_row, reconstructed_row).sum())
        similarity = 1.0 - hamming / float(variant_count)
        precision, recall, row_f1, tp, fp, fn = _row_f1(true_row, reconstructed_row)
        out.append(
            {
                "beacon_size": beacon_size,
                "window_id": window_id,
                "replicate_id": replicate_id,
                "method": method,
                "reference_population": reference_population,
                "template_design": template_design,
                "panel_size": panel_size,
                "beacon_sample_id": beacon_ids[int(true_index)],
                "true_row_index": int(true_index),
                "matched_reconstructed_row_index": int(reconstructed_index),
                "variant_count": variant_count,
                "best_match_hamming": hamming,
                "best_match_similarity": similarity,
                "row_precision": precision,
                "row_recall": recall,
                "row_f1": row_f1,
                "true_carrier_count": int(true_row.sum()),
                "reconstructed_carrier_count": int(reconstructed_row.sum()),
                "true_positives": tp,
                "false_positives": fp,
                "false_negatives": fn,
                "max_count_error": max_count_error,
                "source_overlap_count": source_overlap_count,
                "invalid_run": invalid_run,
            }
        )
    return out


def _aggregate_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    keys = sorted(
        {
            (
                int(row["beacon_size"]),
                row["method"],
                row["reference_population"],
                row["template_design"],
                str(row["panel_size"]),
            )
            for row in rows
        }
    )
    for beacon_size, method, reference_population, template_design, panel_size in keys:
        selected = [
            row
            for row in rows
            if int(row["beacon_size"]) == beacon_size
            and row["method"] == method
            and row["reference_population"] == reference_population
            and row["template_design"] == template_design
            and str(row["panel_size"]) == panel_size
        ]
        similarities = np.asarray([float(row["best_match_similarity"]) for row in selected], dtype=float)
        hamming = np.asarray([float(row["best_match_hamming"]) for row in selected], dtype=float)
        row_f1 = np.asarray([float(row["row_f1"]) for row in selected], dtype=float)
        aggregate = {
            "beacon_size": beacon_size,
            "method": method,
            "reference_population": reference_population,
            "template_design": template_design,
            "panel_size": "" if panel_size == "None" else panel_size,
            "n_individual_rows": len(selected),
            "n_runs": len({(row["window_id"], row["replicate_id"]) for row in selected}),
            "mean_best_match_hamming": float(np.mean(hamming)),
            "median_best_match_hamming": float(np.median(hamming)),
            "p90_best_match_hamming": float(np.quantile(hamming, 0.90)),
            "mean_best_match_similarity": float(np.mean(similarities)),
            "median_best_match_similarity": float(np.median(similarities)),
            "mean_row_f1": float(np.nanmean(row_f1)),
            "median_row_f1": float(np.nanmedian(row_f1)),
            "total_source_overlap_count": int(sum(int(row["source_overlap_count"]) for row in selected)),
            "invalid_individual_rows": int(sum(bool(row["invalid_run"]) for row in selected)),
        }
        for threshold in SIMILARITY_THRESHOLDS:
            label = str(threshold).replace(".", "")
            aggregate[f"proportion_similarity_ge_{label}"] = float(np.mean(similarities >= threshold))
        out.append(aggregate)
    return out


def _write_report(path: Path, aggregate_rows: list[dict[str, Any]], config: dict[str, Any]) -> None:
    lines = [
        "# Genomic Data-Sharing Participant Privacy Metrics",
        "",
        "## Configuration",
        "",
        f"- windows: {config['max_windows']}",
        f"- replicates: {config['max_replicates']}",
        f"- beacon sizes: {','.join(str(size) for size in config['beacon_sizes'])}",
        f"- panel sizes: {','.join(str(size) for size in config['panel_sizes'])}",
        f"- reference panel mode: {config['reference_panel_mode']}",
        f"- pairwise budget: {config['epochs']} epochs x {config['proposals_per_epoch']} proposals",
        "",
        "## Aggregate Participant-Level Risk",
        "",
        "| N | method | n_runs | rows | mean_similarity | median_similarity | mean_hamming | row_f1 | sim>=0.90 | sim>=0.95 | sim>=0.99 | invalid_rows | overlap |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['beacon_size']} | {row['method']} | {row['n_runs']} | {row['n_individual_rows']} | "
            f"{row['mean_best_match_similarity']:.6f} | {row['median_best_match_similarity']:.6f} | "
            f"{row['mean_best_match_hamming']:.2f} | {row['mean_row_f1']:.6f} | "
            f"{row['proportion_similarity_ge_09']:.6f} | {row['proportion_similarity_ge_095']:.6f} | "
            f"{row['proportion_similarity_ge_099']:.6f} | {row['invalid_individual_rows']} | "
            f"{row['total_source_overlap_count']} |"
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "These metrics describe post-hoc participant-level reconstruction risk after row alignment. "
            "They do not add reconstruction-time information; they only reinterpret the same reconstructed "
            "matrices at the individual level.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--beacon-size", type=int, default=25)
    parser.add_argument(
        "--beacon-sizes",
        default=None,
        help="Comma-separated beacon sizes (default: --beacon-size only).",
    )
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--proposals-per-epoch", type=int, default=5000)
    parser.add_argument("--early-stopping", type=int, default=5)
    parser.add_argument("--panel-sizes", default="25,50,100")
    parser.add_argument(
        "--shared-reference-panel",
        action="store_true",
        help="Use one R=100 reference panel for LD/Saleem/pairwise and template rows (disjoint from beacon only).",
    )
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument(
        "--merge-existing",
        action="store_true",
        help="Keep rows in out-dir CSVs that are not replaced by this run.",
    )
    add_saleem_arguments(parser)
    args = parser.parse_args()

    out_dir = args.out_dir or (DEFAULT_SHARED_OUT_DIR if args.shared_reference_panel else OUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    windows = _read_csv(MULTIWINDOW_DIR / "multiwindow_manifest.csv")[: max(0, args.max_windows)]
    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]
    panel_sizes = sorted({int(value) for value in args.panel_sizes.split(",") if value.strip()} | {100})
    beacon_sizes = sorted(
        {int(value) for value in (args.beacon_sizes or str(args.beacon_size)).split(",") if value.strip()}
    )
    config = {
        "max_windows": args.max_windows,
        "max_replicates": args.max_replicates,
        "beacon_sizes": beacon_sizes,
        "epochs": args.epochs,
        "proposals_per_epoch": args.proposals_per_epoch,
        "early_stopping": args.early_stopping,
        "panel_sizes": panel_sizes,
        "reference_panel_mode": "shared" if args.shared_reference_panel else "dual_disjoint",
        **saleem_config_dict(args),
    }

    row_metrics: list[dict[str, Any]] = []
    run_summaries: list[dict[str, Any]] = []
    for beacon_size in beacon_sizes:
        for window in windows:
            window_id = window["window_id"]
            sample_ids, variant_ids, positions, target_alleles, carrier_matrix = load_raw_carriers(
                Path(window["raw_path"])
            )
            wrong_ids, wrong_variant_ids, _wrong_positions, _wrong_targets, wrong_matrix = load_raw_carriers(
                MULTIWINDOW_DIR / "wrong_population" / f"{window_id}.raw",
                target_alleles=target_alleles,
            )
            if wrong_variant_ids != variant_ids:
                raise ValueError(f"Wrong-population variant order mismatch for {window_id}.")

            for replicate_id in replicate_ids:
                seeds = seeds_by_replicate[replicate_id]
                split = split_samples(
                    sample_rows,
                    replicate_id,
                    beacon_size,
                    shared_reference_panel=args.shared_reference_panel,
                )
                template_ref_ids = split["public_ref"] if args.shared_reference_panel else split["arg_ref"]
                beacon_matrix = matrix_for_samples(sample_ids, carrier_matrix, split["beacon"])
                public_ref_matrix = matrix_for_samples(sample_ids, carrier_matrix, split["public_ref"])
                arg_ref_matrix = matrix_for_samples(sample_ids, carrier_matrix, template_ref_ids)
                wrong_ref_ids = wrong_ids[: len(template_ref_ids)]
                wrong_ref_matrix = matrix_for_samples(wrong_ids, wrong_matrix, wrong_ref_ids)
                target_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                reference_pairwise = pairwise_moments(public_ref_matrix)
                split_id = f"{window_id}_replicate_{replicate_id}_N{beacon_size}"
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
                        metadata={"window_id": window_id, "experiment": "participant_privacy"},
                    ),
                )
                truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(split["beacon"]))
                runtime_config = {
                    "window_id": window_id,
                    "replicate_id": replicate_id,
                    "beacon_size": beacon_size,
                    "config": config,
                }
                run_specs: list[dict[str, Any]] = []
                if not args.saleem_only:
                    run_specs.append(
                        {
                            "method": "independent_count",
                            "reference_population": "none",
                            "template_design": "count_only",
                            "panel_size": None,
                            "reference_ids": [],
                            "result": _with_config(
                                independent_count(problem, reference_pairwise=reference_pairwise),
                                config=runtime_config,
                                optimizer_seed=seeds["optimizer_seed"],
                            ),
                        }
                    )
                    run_specs.append(
                        {
                            "method": "pairwise_refinement",
                            "reference_population": "eur",
                            "template_design": "pairwise_ld",
                            "panel_size": len(split["public_ref"]),
                            "reference_ids": split["public_ref"],
                            "result": _with_config(
                                pairwise_refinement(
                                    problem,
                                    reference_pairwise,
                                    epochs=args.epochs,
                                    proposals_per_epoch=args.proposals_per_epoch,
                                    early_stopping=args.early_stopping,
                                ),
                                config=runtime_config,
                                optimizer_seed=seeds["optimizer_seed"],
                            ),
                        }
                    )
                if args.run_saleem:
                    run_specs.append(
                        {
                            "method": SALEEM_RECONSTRUCTION_ATTACK,
                            "reference_population": "eur",
                            "template_design": "saleem_ld",
                            "panel_size": len(split["public_ref"]),
                            "reference_ids": split["public_ref"],
                            "result": _with_config(
                                run_saleem_attack(problem, public_ref_matrix, beacon_matrix, args),
                                config=runtime_config,
                                optimizer_seed=seeds["optimizer_seed"],
                            ),
                        }
                    )

                if not args.saleem_only:
                    for panel_size in panel_sizes:
                        panel_ids = sample_reference_ids(
                            template_ref_ids,
                            panel_size=panel_size,
                            seed=seeds["group_selection_seed"] + panel_size,
                        )
                        panel_matrix = matrix_for_samples(sample_ids, carrier_matrix, panel_ids)
                        run_specs.append(
                            {
                                "method": f"matched_template_n{panel_size}",
                                "reference_population": "eur",
                                "template_design": "matched",
                                "panel_size": panel_size,
                                "reference_ids": panel_ids,
                                "result": _with_config(
                                    reference_haplotype_reconstruction(problem, panel_matrix, metadata=runtime_config),
                                    config=runtime_config,
                                    optimizer_seed=seeds["optimizer_seed"],
                                ),
                            }
                        )

                    shuffled_matrix = shuffle_reference_columns(
                        arg_ref_matrix,
                        seed=seeds["group_selection_seed"] + 19_001,
                    )
                    run_specs.extend(
                        [
                            {
                                "method": "shuffled_template_n100",
                                "reference_population": "eur",
                                "template_design": "column_shuffled",
                                "panel_size": len(template_ref_ids),
                                "reference_ids": template_ref_ids,
                                "result": _with_config(
                                    reference_haplotype_reconstruction(problem, shuffled_matrix, metadata=runtime_config),
                                    config=runtime_config,
                                    optimizer_seed=seeds["optimizer_seed"],
                                ),
                            },
                            {
                                "method": "wrong_population_template_n100",
                                "reference_population": "afr",
                                "template_design": "wrong_population",
                                "panel_size": len(wrong_ref_ids),
                                "reference_ids": wrong_ref_ids,
                                "result": _with_config(
                                    reference_haplotype_reconstruction(problem, wrong_ref_matrix, metadata=runtime_config),
                                    config=runtime_config,
                                    optimizer_seed=seeds["optimizer_seed"],
                                ),
                            },
                        ]
                    )

                for spec in run_specs:
                    source_overlap_count = len(set(split["beacon"]).intersection(spec["reference_ids"]))
                    evaluation = evaluate_reconstruction(
                        spec["result"],
                        truth,
                        target_counts=target_counts,
                        reference_pairwise=reference_pairwise,
                        source_sets={
                            "beacon": split["beacon"],
                            "template_reference": spec["reference_ids"],
                        },
                        accepted_manifest_ids=None,
                    )
                    count_metrics = compute_count_metrics(spec["result"].reconstructed_binary, target_counts)
                    rows = _row_privacy_rows(
                        beacon_size=beacon_size,
                        window_id=window_id,
                        replicate_id=replicate_id,
                        method=spec["method"],
                        reference_population=spec["reference_population"],
                        template_design=spec["template_design"],
                        panel_size=spec["panel_size"],
                        beacon_ids=split["beacon"],
                        true_matrix=beacon_matrix,
                        reconstructed_matrix=spec["result"].reconstructed_binary,
                        max_count_error=count_metrics.max_count_error,
                        source_overlap_count=source_overlap_count,
                        invalid_run=evaluation.invalid_run.is_invalid,
                    )
                    row_metrics.extend(rows)
                    run_summaries.append(
                        {
                            "beacon_size": beacon_size,
                            "window_id": window_id,
                            "replicate_id": replicate_id,
                            "method": spec["method"],
                            "f1": evaluation.metrics.f1,
                            "alignment_total_hamming_cost": evaluation.metrics.alignment_total_hamming_cost,
                            "max_count_error": count_metrics.max_count_error,
                            "source_overlap_count": source_overlap_count,
                            "invalid_run": evaluation.invalid_run.is_invalid,
                            "invalid_reasons": ";".join(evaluation.invalid_run.reasons),
                        }
                    )
                print(f"completed {split_id}")

    run_summaries = merge_saleem_rows(
        new_rows=run_summaries,
        merge_metrics=args.merge_metrics,
        saleem_only=args.saleem_only,
    )
    row_metrics = merge_saleem_rows(
        new_rows=row_metrics,
        merge_metrics=out_dir / "row_privacy_metrics.csv",
        saleem_only=args.saleem_only,
    )
    if args.merge_existing:
        run_summaries = _merge_existing_rows(
            run_summaries,
            metrics_path=out_dir / "row_privacy_run_summary.csv",
            key_fn=_run_summary_key,
        )
        row_metrics = _merge_existing_rows(
            row_metrics,
            metrics_path=out_dir / "row_privacy_metrics.csv",
            key_fn=_row_metric_key,
        )
    aggregate = _aggregate_rows(row_metrics)
    _write_csv(out_dir / "row_privacy_metrics.csv", row_metrics)
    _write_csv(out_dir / "row_privacy_aggregate.csv", aggregate)
    _write_csv(out_dir / "row_privacy_run_summary.csv", run_summaries)
    _write_json(out_dir / "row_privacy_manifest.json", {"manifest_id": MANIFEST_ID, "config": config})
    _write_report(out_dir / "row_privacy_report.md", aggregate, config)
    print(f"wrote {out_dir / 'row_privacy_report.md'}")


if __name__ == "__main__":
    main()
