"""Run cohort-size sensitivity for Beacon reconstruction experiments."""

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
from beacon_reconstruction.ld_af_baseline import SALEEM_RECONSTRUCTION_ATTACK
from beacon_reconstruction.core import ExperimentConfig, BeaconAttackProblem, BeaconEvaluationTruth, ReconstructionResult
from beacon_reconstruction.evaluation import EvaluationResult, evaluate_reconstruction
from scripts.ld_af_helpers import (
    add_saleem_arguments,
    merge_saleem_overlap_rows,
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
from scripts.experiment_runtime_helpers import (
    attach_reference_panel_mode,
    run_timed,
    runtime_metadata_fields,
)


MANIFEST_ID = "1000g_eur_hg38_cohort_scaling_v1"
BASELINE_METHOD = SALEEM_RECONSTRUCTION_ATTACK
DEFAULT_MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
DEFAULT_OUT_DIR = ROOT / "results" / "experiments" / "n_sensitivity"
DEFAULT_SHARED_OUT_DIR = ROOT / "results" / "experiments" / "n_sensitivity_shared_panel"


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


def _record(
    *,
    beacon_size: int,
    window_id: str,
    replicate_id: str,
    method: str,
    reference_population: str,
    template_design: str,
    panel_size: int | None,
    result: ReconstructionResult,
    evaluation: EvaluationResult,
    source_overlap_count: int,
    runtime_fields: dict[str, Any] | None = None,
) -> dict[str, Any]:
    metrics = evaluation.metrics
    row = {
        "beacon_size": beacon_size,
        "window_id": window_id,
        "replicate_id": replicate_id,
        "method": method,
        "method_name": result.method_name,
        "reference_population": reference_population,
        "template_design": template_design,
        "panel_size": panel_size,
        "source_overlap_count": source_overlap_count,
        "precision": metrics.precision,
        "recall": metrics.recall,
        "f1": metrics.f1,
        "pairwise_error": metrics.pairwise_error,
        "alignment_total_hamming_cost": metrics.alignment_total_hamming_cost,
        "max_count_error": metrics.max_count_error,
        "invalid_run": evaluation.invalid_run.is_invalid,
        "invalid_reasons": ";".join(evaluation.invalid_run.reasons),
    }
    if runtime_fields:
        row.update(runtime_fields)
    return row


def _overlap_row(
    *,
    beacon_size: int,
    window_id: str,
    replicate_id: str,
    method: str,
    beacon_ids: list[str],
    reference_ids: list[str],
) -> dict[str, Any]:
    overlap = sorted(set(beacon_ids).intersection(reference_ids))
    return {
        "beacon_size": beacon_size,
        "window_id": window_id,
        "replicate_id": replicate_id,
        "method": method,
        "beacon_sample_count": len(beacon_ids),
        "reference_sample_count": len(reference_ids),
        "overlap_count": len(overlap),
        "overlap_sample_ids": ";".join(overlap),
    }


def _aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    keys = sorted({(int(row["beacon_size"]), row["method"]) for row in rows})
    for beacon_size, method in keys:
        selected = [row for row in rows if int(row["beacon_size"]) == beacon_size and row["method"] == method]
        out.append(
            {
                "beacon_size": beacon_size,
                "method": method,
                "reference_population": selected[0]["reference_population"],
                "template_design": selected[0]["template_design"],
                "panel_size": selected[0]["panel_size"],
                "n_runs": len(selected),
                "mean_f1": float(np.mean([float(row["f1"]) for row in selected])),
                "median_f1": float(np.median([float(row["f1"]) for row in selected])),
                "mean_hamming": float(np.mean([float(row["alignment_total_hamming_cost"]) for row in selected])),
                "mean_pairwise_error": float(np.mean([float(row["pairwise_error"]) for row in selected])),
                "total_source_overlap_count": int(sum(int(row["source_overlap_count"]) for row in selected)),
                "invalid_run_count": int(sum(bool(row["invalid_run"]) for row in selected)),
            }
        )
        if "wall_time_sec" in selected[0]:
            wall_times = [float(row["wall_time_sec"]) for row in selected]
            out[-1]["mean_wall_time_sec"] = float(np.mean(wall_times))
            out[-1]["median_wall_time_sec"] = float(np.median(wall_times))
        if "saleem_outer_iterations" in selected[0]:
            saleem_iters = [
                float(row["saleem_outer_iterations"])
                for row in selected
                if row.get("saleem_outer_iterations") not in ("", None)
            ]
            if saleem_iters:
                out[-1]["mean_saleem_outer_iterations"] = float(np.mean(saleem_iters))
    return out


def _pairwise_deltas(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {
        (int(row["beacon_size"]), row["window_id"], row["replicate_id"], row["method"]): row for row in rows
    }
    units = sorted({(int(row["beacon_size"]), row["window_id"], row["replicate_id"]) for row in rows})
    methods = sorted({row["method"] for row in rows} - {BASELINE_METHOD})
    out: list[dict[str, Any]] = []
    for beacon_size, window_id, replicate_id in units:
        baseline = by_key[(beacon_size, window_id, replicate_id, BASELINE_METHOD)]
        for method in methods:
            row = by_key[(beacon_size, window_id, replicate_id, method)]
            out.append(
                {
                    "beacon_size": beacon_size,
                    "window_id": window_id,
                    "replicate_id": replicate_id,
                    "method": method,
                    "delta_f1_vs_saleem": float(row["f1"]) - float(baseline["f1"]),
                    "delta_hamming_vs_saleem": float(row["alignment_total_hamming_cost"])
                    - float(baseline["alignment_total_hamming_cost"]),
                }
            )
    return out


def _write_report(
    path: Path,
    *,
    aggregate_rows: list[dict[str, Any]],
    delta_rows: list[dict[str, Any]],
    overlap_rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> None:
    lines = [
        "# Beacon-Size Sensitivity For Genomic Data-Sharing Privacy Risk",
        "",
        "## Configuration",
        "",
        f"- windows: {config['max_windows']}",
        f"- replicates: {config['max_replicates']}",
        f"- beacon sizes: {','.join(str(size) for size in config['beacon_sizes'])}",
        f"- panel sizes: {','.join(str(size) for size in config['panel_sizes'])}",
        f"- reference panel mode: {config.get('reference_panel_mode', 'dual_disjoint')}",
        f"- Saleem baseline: corr={config['saleem_corr_epoch']}, freq={config['saleem_freq_epoch']}, "
        f"max_iter={config['saleem_max_iterations']}, flip<{config['saleem_flip_threshold']}, "
        f"lr={config['saleem_learning_rate']}",
        "",
        "## Aggregate Metrics",
        "",
        "| N | method | n_runs | mean_f1 | median_f1 | mean_hamming | mean_wall_s | overlap | invalid |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        wall = row.get("mean_wall_time_sec")
        wall_text = f"{wall:.3f}" if wall is not None else "—"
        lines.append(
            f"| {row['beacon_size']} | {row['method']} | {row['n_runs']} | "
            f"{row['mean_f1']:.6f} | {row['median_f1']:.6f} | {row['mean_hamming']:.2f} | "
            f"{wall_text} | {row['total_source_overlap_count']} | {row['invalid_run_count']} |"
        )

    lines.extend(
        [
            "",
            "## Deltas Versus Saleem Baseline",
            "",
            "| N | method | mean_delta_f1 | median_delta_f1 | wins | ties | losses |",
            "|---:|---|---:|---:|---:|---:|---:|",
        ]
    )
    for beacon_size, method in sorted({(int(row["beacon_size"]), row["method"]) for row in delta_rows}):
        selected = [
            row for row in delta_rows if int(row["beacon_size"]) == beacon_size and row["method"] == method
        ]
        deltas = np.asarray([float(row["delta_f1_vs_saleem"]) for row in selected], dtype=float)
        lines.append(
            f"| {beacon_size} | {method} | {float(np.mean(deltas)):.6f} | {float(np.median(deltas)):.6f} | "
            f"{int(np.sum(deltas > 0))} | {int(np.sum(deltas == 0))} | {int(np.sum(deltas < 0))} |"
        )

    lines.extend(
        [
            "",
            "## Source Overlap",
            "",
            f"Total source overlap across all rows: {sum(int(row['overlap_count']) for row in overlap_rows)}.",
            "",
            "## Claim Boundary",
            "",
            "This analysis evaluates whether cohort size changes participant-level genomic reconstruction risk. "
            "The fixed sample split supports `N=25`, `N=50`, `N=100`, and `N=250` with disjoint beacon and "
            "reference groups within each replicate. "
            f"Reference panel mode: `{config.get('reference_panel_mode', 'dual_disjoint')}`.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--beacon-sizes", default="25,50,100,250")
    parser.add_argument("--panel-sizes", default="25,50,100")
    parser.add_argument("--pairwise-epochs", type=int, default=20)
    parser.add_argument("--pairwise-proposals-per-epoch", type=int, default=5000)
    parser.add_argument("--pairwise-early-stopping", type=int, default=5)
    parser.add_argument(
        "--shared-reference-panel",
        action="store_true",
        help="Use one R=100 reference panel for LD/Saleem/pairwise and template rows (disjoint from beacon only).",
    )
    add_saleem_arguments(parser)
    parser.add_argument("--multiwindow-dir", type=Path, default=DEFAULT_MULTIWINDOW_DIR)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()

    out_dir = args.out_dir or (DEFAULT_SHARED_OUT_DIR if args.shared_reference_panel else DEFAULT_OUT_DIR)
    multiwindow_dir = args.multiwindow_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    windows = _read_csv(multiwindow_dir / "multiwindow_manifest.csv")[: max(0, args.max_windows)]
    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]
    beacon_sizes = sorted({int(value) for value in args.beacon_sizes.split(",") if value.strip()})
    panel_sizes = sorted({int(value) for value in args.panel_sizes.split(",") if value.strip()} | {100})
    config = {
        "max_windows": args.max_windows,
        "max_replicates": args.max_replicates,
        "beacon_sizes": beacon_sizes,
        "panel_sizes": panel_sizes,
        "pairwise_epochs": args.pairwise_epochs,
        "pairwise_proposals_per_epoch": args.pairwise_proposals_per_epoch,
        "pairwise_early_stopping": args.pairwise_early_stopping,
        **saleem_config_dict(args),
        "reference_panel_mode": "shared" if args.shared_reference_panel else "dual_disjoint",
        "multiwindow_dir": str(multiwindow_dir),
        "out_dir": str(out_dir),
    }

    rows: list[dict[str, Any]] = []
    overlap_rows: list[dict[str, Any]] = []
    for window in windows:
        window_id = window["window_id"]
        sample_ids, variant_ids, positions, target_alleles, carrier_matrix = load_raw_carriers(
            Path(window["raw_path"])
        )
        wrong_ids, wrong_variant_ids, _wrong_positions, _wrong_targets, wrong_matrix = load_raw_carriers(
            multiwindow_dir / "wrong_population" / f"{window_id}.raw",
            target_alleles=target_alleles,
        )
        if wrong_variant_ids != variant_ids:
            raise ValueError(f"Wrong-population variant order mismatch for {window_id}.")

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
                beacon_matrix = matrix_for_samples(sample_ids, carrier_matrix, split["beacon"])
                public_ref_matrix = matrix_for_samples(sample_ids, carrier_matrix, split["public_ref"])
                arg_ref_matrix = matrix_for_samples(sample_ids, carrier_matrix, template_ref_ids)
                wrong_ref_ids = wrong_ids[: len(template_ref_ids)]
                wrong_ref_matrix = matrix_for_samples(wrong_ids, wrong_matrix, wrong_ref_ids)
                target_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                reference_pairwise = pairwise_moments(public_ref_matrix)
                split_id = f"{window_id}_replicate_{replicate_id}_N{beacon_size}"
                runtime_config = {
                    "window_id": window_id,
                    "replicate_id": replicate_id,
                    "beacon_size": beacon_size,
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
                        metadata={"window_id": window_id, "experiment": "n_sensitivity"},
                    ),
                )
                truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(split["beacon"]))

                run_spec_defs: list[dict[str, Any]] = []
                if not args.saleem_only:
                    run_spec_defs.append(
                        {
                            "method": "independent_count",
                            "reference_population": "none",
                            "template_design": "count_only",
                            "panel_size": None,
                            "reference_ids": [],
                            "run": lambda problem=problem, reference_pairwise=reference_pairwise, runtime_config=runtime_config, seeds=seeds: _with_config(
                                independent_count(problem, reference_pairwise=reference_pairwise),
                                config=runtime_config,
                                optimizer_seed=seeds["optimizer_seed"],
                            ),
                        }
                    )
                    run_spec_defs.append(
                        {
                            "method": "pairwise_refinement",
                            "reference_population": "eur",
                            "template_design": "pairwise_ld",
                            "panel_size": len(split["public_ref"]),
                            "reference_ids": split["public_ref"],
                            "run": lambda problem=problem, reference_pairwise=reference_pairwise, runtime_config=runtime_config, seeds=seeds: _with_config(
                                pairwise_refinement(
                                    problem,
                                    reference_pairwise,
                                    epochs=args.pairwise_epochs,
                                    proposals_per_epoch=args.pairwise_proposals_per_epoch,
                                    early_stopping=args.pairwise_early_stopping,
                                ),
                                config=runtime_config,
                                optimizer_seed=seeds["optimizer_seed"],
                            ),
                        }
                    )
                if args.run_saleem:
                    run_spec_defs.append(
                        {
                            "method": BASELINE_METHOD,
                            "reference_population": "eur",
                            "template_design": "saleem_ld",
                            "panel_size": len(split["public_ref"]),
                            "reference_ids": split["public_ref"],
                            "run": lambda problem=problem, public_ref_matrix=public_ref_matrix, beacon_matrix=beacon_matrix, runtime_config=runtime_config, seeds=seeds: _with_config(
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
                        run_spec_defs.append(
                            {
                                "method": f"matched_template_n{panel_size}",
                                "reference_population": "eur",
                                "template_design": "matched",
                                "panel_size": panel_size,
                                "reference_ids": panel_ids,
                                "run": lambda problem=problem, panel_matrix=panel_matrix, runtime_config=runtime_config, seeds=seeds: _with_config(
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
                    run_spec_defs.extend(
                        [
                            {
                                "method": "shuffled_template_n100",
                                "reference_population": "eur",
                                "template_design": "column_shuffled",
                                "panel_size": len(template_ref_ids),
                                "reference_ids": template_ref_ids,
                                "run": lambda problem=problem, shuffled_matrix=shuffled_matrix, runtime_config=runtime_config, seeds=seeds: _with_config(
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
                                "run": lambda problem=problem, wrong_ref_matrix=wrong_ref_matrix, runtime_config=runtime_config, seeds=seeds: _with_config(
                                    reference_haplotype_reconstruction(problem, wrong_ref_matrix, metadata=runtime_config),
                                    config=runtime_config,
                                    optimizer_seed=seeds["optimizer_seed"],
                                ),
                            },
                        ]
                    )

                for spec_def in run_spec_defs:
                    result, wall_time_sec = run_timed(spec_def["run"])
                    result = attach_reference_panel_mode(result, shared=args.shared_reference_panel)
                    runtime_fields = runtime_metadata_fields(result, wall_time_sec)
                    source_overlap_count = len(set(split["beacon"]).intersection(spec_def["reference_ids"]))
                    overlap_rows.append(
                        _overlap_row(
                            beacon_size=beacon_size,
                            window_id=window_id,
                            replicate_id=replicate_id,
                            method=spec_def["method"],
                            beacon_ids=split["beacon"],
                            reference_ids=spec_def["reference_ids"],
                        )
                    )
                    evaluation = evaluate_reconstruction(
                        result,
                        truth,
                        target_counts=target_counts,
                        reference_pairwise=reference_pairwise,
                        source_sets={
                            "beacon": split["beacon"],
                            "template_reference": spec_def["reference_ids"],
                        },
                        accepted_manifest_ids=None,
                    )
                    rows.append(
                        _record(
                            beacon_size=beacon_size,
                            window_id=window_id,
                            replicate_id=replicate_id,
                            method=spec_def["method"],
                            reference_population=spec_def["reference_population"],
                            template_design=spec_def["template_design"],
                            panel_size=spec_def["panel_size"],
                            result=result,
                            evaluation=evaluation,
                            source_overlap_count=source_overlap_count,
                            runtime_fields=runtime_fields,
                        )
                    )
                print(f"completed {split_id}")

    rows = merge_saleem_rows(new_rows=rows, merge_metrics=args.merge_metrics, saleem_only=args.saleem_only)
    overlap_rows = merge_saleem_overlap_rows(
        new_rows=overlap_rows,
        merge_metrics=args.merge_metrics,
        saleem_only=args.saleem_only,
        overlap_filename="n_sensitivity_source_overlap.csv",
    )
    aggregate = _aggregate(rows)
    deltas = _pairwise_deltas(rows)
    _write_csv(out_dir / "n_sensitivity_metrics.csv", rows)
    _write_csv(out_dir / "n_sensitivity_aggregate.csv", aggregate)
    _write_csv(out_dir / "n_sensitivity_pairwise_deltas.csv", deltas)
    _write_csv(out_dir / "n_sensitivity_source_overlap.csv", overlap_rows)
    _write_json(out_dir / "n_sensitivity_manifest.json", {"manifest_id": MANIFEST_ID, "config": config})
    _write_report(
        out_dir / "n_sensitivity_report.md",
        aggregate_rows=aggregate,
        delta_rows=deltas,
        overlap_rows=overlap_rows,
        config=config,
    )
    print(f"wrote {out_dir / 'n_sensitivity_report.md'}")


if __name__ == "__main__":
    main()
