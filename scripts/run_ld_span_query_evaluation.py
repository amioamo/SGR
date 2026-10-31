"""Evaluate reconstruction across LD, density, and genomic-span query regimes.

Each query set has the same SNP count. Regimes differ in residual linkage,
inter-SNP spacing, or genomic span, and sit in additional regions outside the
six main dense windows. The default design uses one shared EUR reference panel.
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

from beacon_reconstruction.baselines import independent_count, pairwise_moments, pairwise_refinement
from beacon_reconstruction.reconstruction import TEMPLATE_SELECTION_RULE, reference_haplotype_reconstruction
from beacon_reconstruction.core import ExperimentConfig, BeaconAttackProblem, BeaconEvaluationTruth, ReconstructionResult
from beacon_reconstruction.evaluation import EvaluationResult, evaluate_reconstruction
from beacon_reconstruction.ld_af_baseline import SALEEM_RECONSTRUCTION_ATTACK
from scripts.experiment_runtime_helpers import attach_reference_panel_mode, run_timed, runtime_metadata_fields
from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers, matrix_for_samples as _matrix_for_samples
from scripts.lib.sample_splits import load_sample_manifest, load_seed_manifest, split_samples
from scripts.lib.reference_helpers import shuffle_reference_columns as _shuffle_template_columns
from scripts.ld_af_helpers import (
    add_saleem_arguments,
    merge_saleem_overlap_rows,
    merge_saleem_rows,
    run_saleem_attack,
    saleem_config_dict,
)


MANIFEST_ID = "1000g_eur_hg38_ld_span_queries_v1"
DEFAULT_MANIFEST = ROOT / "data" / "ld_span_queries" / "ld_span_query_manifest.csv"
DEFAULT_OUT_DIR = ROOT / "results" / "block_inheritance_v1" / "ld_span_queries_shared_panel"
OVERLAP_FILENAME = "ld_span_query_source_overlap.csv"


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
    for row in rows:
        for key in row:
            if key not in fieldnames:
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
    query_set: dict[str, str],
    beacon_size: int,
    replicate_id: str,
    method: str,
    reference_population: str,
    template_design: str,
    result: ReconstructionResult,
    evaluation: EvaluationResult,
    source_overlap_count: int,
    runtime_fields: dict[str, Any],
) -> dict[str, Any]:
    metrics = evaluation.metrics
    row = {
        "query_set_id": query_set["query_set_id"],
        "regime": query_set["regime"],
        "chrom": query_set["chrom"],
        "region_index": query_set["region_index"],
        "variant_count": int(query_set["variant_count"]),
        "span_bp": int(float(query_set["span_bp"])),
        "median_inter_snp_bp": float(query_set["median_inter_snp_bp"]),
        "snps_per_mb": float(query_set["snps_per_mb"]),
        "mean_abs_corr": float(query_set["mean_abs_corr"]),
        "max_pairwise_r2": float(query_set["max_pairwise_r2"]),
        "beacon_size": beacon_size,
        "replicate_id": replicate_id,
        "method": method,
        "method_name": result.method_name,
        "reference_population": reference_population,
        "template_design": template_design,
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
    row.update(runtime_fields)
    return row


def _aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = sorted({(int(row["beacon_size"]), row["regime"], row["method"]) for row in rows})
    out: list[dict[str, Any]] = []
    for beacon_size, regime, method in keys:
        selected = [
            row
            for row in rows
            if int(row["beacon_size"]) == beacon_size and row["regime"] == regime and row["method"] == method
        ]
        out.append(
            {
                "beacon_size": beacon_size,
                "regime": regime,
                "method": method,
                "n_runs": len(selected),
                "mean_span_mb": float(np.mean([float(row["span_bp"]) for row in selected])) / 1e6,
                "mean_median_gap_kb": float(np.mean([float(row["median_inter_snp_bp"]) for row in selected])) / 1e3,
                "mean_snps_per_mb": float(np.mean([float(row["snps_per_mb"]) for row in selected])),
                "mean_abs_corr": float(np.mean([float(row["mean_abs_corr"]) for row in selected])),
                "mean_max_pairwise_r2": float(np.mean([float(row["max_pairwise_r2"]) for row in selected])),
                "mean_f1": float(np.mean([float(row["f1"]) for row in selected])),
                "median_f1": float(np.median([float(row["f1"]) for row in selected])),
                "total_source_overlap_count": int(sum(int(row["source_overlap_count"]) for row in selected)),
                "invalid_run_count": int(sum(bool(row["invalid_run"]) for row in selected)),
            }
        )
    return out


def _pairwise_deltas(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {
        (row["query_set_id"], str(row["replicate_id"]), int(row["beacon_size"]), row["method"]): row for row in rows
    }
    units = sorted({(row["query_set_id"], str(row["replicate_id"]), int(row["beacon_size"])) for row in rows})
    methods = sorted({row["method"] for row in rows} - {"pairwise_refinement"})
    out: list[dict[str, Any]] = []
    for query_set_id, replicate_id, beacon_size in units:
        pairwise = by_key.get((query_set_id, replicate_id, beacon_size, "pairwise_refinement"))
        if pairwise is None:
            continue
        for method in methods:
            row = by_key.get((query_set_id, replicate_id, beacon_size, method))
            if row is None:
                continue
            out.append(
                {
                    "query_set_id": query_set_id,
                    "regime": row["regime"],
                    "chrom": row["chrom"],
                    "region_index": row["region_index"],
                    "beacon_size": beacon_size,
                    "replicate_id": replicate_id,
                    "method": method,
                    "delta_f1_vs_pairwise": float(row["f1"]) - float(pairwise["f1"]),
                }
            )
    return out


def _dense_deltas(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    dense = {
        (row["chrom"], str(row["region_index"]), str(row["replicate_id"]), int(row["beacon_size"]), row["method"]): row
        for row in rows
        if row["regime"] == "dense"
    }
    out: list[dict[str, Any]] = []
    for row in rows:
        if row["regime"] == "dense":
            continue
        baseline = dense.get(
            (row["chrom"], str(row["region_index"]), str(row["replicate_id"]), int(row["beacon_size"]), row["method"])
        )
        if baseline is None:
            continue
        out.append(
            {
                "query_set_id": row["query_set_id"],
                "regime": row["regime"],
                "chrom": row["chrom"],
                "region_index": row["region_index"],
                "beacon_size": int(row["beacon_size"]),
                "replicate_id": row["replicate_id"],
                "method": row["method"],
                "delta_f1_vs_dense": float(row["f1"]) - float(baseline["f1"]),
            }
        )
    return out


def _fmt(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def _write_report(
    path: Path,
    *,
    aggregate_rows: list[dict[str, Any]],
    pairwise_rows: list[dict[str, Any]],
    dense_rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> None:
    lines = [
        "# LD, Density, And Genomic-Span Query Evaluation",
        "",
        "Fixed SNP count across additional genomic regions. `dense` is the short-span, high-density reference. "
        "`ld_r2_080`, `ld_r2_050`, and `ld_r2_010` cap pairwise carrier r^2 at 0.80, 0.50, and 0.10. `sparse` enforces a 50 kb minimum gap. "
        "`long_range` spaces the same SNP count across 10 Mb.",
        "",
        "## Configuration",
        "",
        f"- query sets: {config['query_set_count']}",
        f"- replicates: {config['max_replicates']}",
        f"- beacon sizes: {config['beacon_sizes']}",
        f"- reference panel: {config['reference_panel_mode']}",
        f"- template selection: {config['template_selection_rule']}",
        "",
        "## Realized Panel Structure And Carrier-Level F1",
        "",
        "| N | regime | method | n | span Mb | gap kb | SNPs/Mb | mean |r| | max r2 | mean F1 | invalid |",
        "|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['beacon_size']} | {row['regime']} | {row['method']} | {row['n_runs']} | "
            f"{_fmt(row['mean_span_mb'])} | {_fmt(row['mean_median_gap_kb'])} | {_fmt(row['mean_snps_per_mb'], 1)} | "
            f"{_fmt(row['mean_abs_corr'])} | {_fmt(row['mean_max_pairwise_r2'])} | {_fmt(row['mean_f1'], 3)} | "
            f"{row['invalid_run_count']} |"
        )

    lines.extend(
        [
            "",
            "## Delta F1 Versus Local LD Refinement",
            "",
            "| N | regime | method | mean delta | median delta | wins | ties | losses |",
            "|---:|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    groups = sorted({(int(row["beacon_size"]), row["regime"], row["method"]) for row in pairwise_rows})
    for beacon_size, regime, method in groups:
        selected = [
            row
            for row in pairwise_rows
            if int(row["beacon_size"]) == beacon_size and row["regime"] == regime and row["method"] == method
        ]
        deltas = np.asarray([float(row["delta_f1_vs_pairwise"]) for row in selected], dtype=float)
        lines.append(
            f"| {beacon_size} | {regime} | {method} | {float(np.mean(deltas)):.3f} | "
            f"{float(np.median(deltas)):.3f} | {int(np.sum(deltas > 0))} | {int(np.sum(np.isclose(deltas, 0)))} | "
            f"{int(np.sum(deltas < 0))} |"
        )

    lines.extend(
        [
            "",
            "## Delta F1 Versus The Dense Panel In The Same Additional Region",
            "",
            "| N | regime | method | mean delta | median delta | wins | ties | losses |",
            "|---:|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    dense_groups = sorted({(int(row["beacon_size"]), row["regime"], row["method"]) for row in dense_rows})
    for beacon_size, regime, method in dense_groups:
        selected = [
            row
            for row in dense_rows
            if int(row["beacon_size"]) == beacon_size and row["regime"] == regime and row["method"] == method
        ]
        deltas = np.asarray([float(row["delta_f1_vs_dense"]) for row in selected], dtype=float)
        lines.append(
            f"| {beacon_size} | {regime} | {method} | {float(np.mean(deltas)):.3f} | "
            f"{float(np.median(deltas)):.3f} | {int(np.sum(deltas > 0))} | {int(np.sum(np.isclose(deltas, 0)))} | "
            f"{int(np.sum(deltas < 0))} |"
        )

    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "These panels are additional regions on chromosomes 1, 20, and 22, disjoint from the six main dense windows. "
            "They vary query structure under a matched shared reference panel. They are not genome-wide random queries, "
            "empirical Beacon logs, or a replacement for the main dense-window analysis.",
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--max-query-sets", type=int, default=0, help="0 means all query sets.")
    parser.add_argument("--regimes", default="", help="Comma-separated regime filter. Empty means all.")
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--beacon-sizes", default="25,50")
    parser.add_argument("--pairwise-epochs", type=int, default=20)
    parser.add_argument("--pairwise-proposals-per-epoch", type=int, default=5000)
    parser.add_argument("--pairwise-early-stopping", type=int, default=5)
    parser.add_argument("--shared-reference-panel", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    add_saleem_arguments(parser)
    args = parser.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    query_sets = _read_csv(args.manifest)
    if args.regimes.strip():
        wanted = {item.strip() for item in args.regimes.split(",") if item.strip()}
        query_sets = [row for row in query_sets if row["regime"] in wanted]
    if args.max_query_sets > 0:
        query_sets = query_sets[: args.max_query_sets]
    if not query_sets:
        raise ValueError("No query sets selected.")

    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]
    beacon_sizes = [int(value) for value in args.beacon_sizes.split(",") if value.strip()]
    config = {
        "manifest": str(args.manifest),
        "query_set_count": len(query_sets),
        "max_replicates": args.max_replicates,
        "beacon_sizes": beacon_sizes,
        "pairwise_epochs": args.pairwise_epochs,
        "pairwise_proposals_per_epoch": args.pairwise_proposals_per_epoch,
        "pairwise_early_stopping": args.pairwise_early_stopping,
        "reference_panel_mode": "shared" if args.shared_reference_panel else "dual_disjoint",
        "template_selection_rule": TEMPLATE_SELECTION_RULE,
        **saleem_config_dict(args),
        "out_dir": str(out_dir),
    }

    rows: list[dict[str, Any]] = []
    overlap_rows: list[dict[str, Any]] = []
    for query_set in query_sets:
        sample_ids, variant_ids, positions, _target_alleles, carrier_matrix = _load_raw_as_target_carriers(
            Path(query_set["raw_path"])
        )
        if len(variant_ids) != int(query_set["variant_count"]):
            raise ValueError(f"Variant count mismatch for {query_set['query_set_id']}.")
        for beacon_size in beacon_sizes:
            for replicate_id in replicate_ids:
                seeds = seeds_by_replicate[replicate_id]
                split = split_samples(
                    sample_rows,
                    replicate_id,
                    beacon_size,
                    shared_reference_panel=args.shared_reference_panel,
                )
                beacon_matrix = _matrix_for_samples(sample_ids, carrier_matrix, split["beacon"])
                public_ref_matrix = _matrix_for_samples(sample_ids, carrier_matrix, split["public_ref"])
                target_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                reference_pairwise = pairwise_moments(public_ref_matrix)
                split_id = f"{query_set['query_set_id']}_replicate_{replicate_id}_N{beacon_size}"
                runtime_config = {
                    "query_set_id": query_set["query_set_id"],
                    "regime": query_set["regime"],
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
                        metadata={"query_set_id": query_set["query_set_id"], "experiment": "ld_span_queries"},
                    ),
                )
                truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(split["beacon"]))
                shuffled_matrix = _shuffle_template_columns(
                    public_ref_matrix,
                    seed=seeds["group_selection_seed"] + 19_001,
                )
                run_specs: list[dict[str, Any]] = []
                if not args.saleem_only:
                    run_specs.extend(
                        [
                            {
                                "method": "independent_count",
                                "reference_population": "none",
                                "template_design": "count_only",
                                "reference_ids": [],
                                "run": lambda problem=problem, reference_pairwise=reference_pairwise, runtime_config=runtime_config, seeds=seeds: _with_config(
                                    independent_count(problem, reference_pairwise=reference_pairwise),
                                    config=runtime_config,
                                    optimizer_seed=seeds["optimizer_seed"],
                                ),
                            },
                            {
                                "method": "pairwise_refinement",
                                "reference_population": "eur",
                                "template_design": "pairwise_ld",
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
                            },
                            {
                                "method": "matched_template_n100",
                                "reference_population": "eur",
                                "template_design": "matched",
                                "reference_ids": split["public_ref"],
                                "run": lambda problem=problem, public_ref_matrix=public_ref_matrix, runtime_config=runtime_config, seeds=seeds: _with_config(
                                    reference_haplotype_reconstruction(problem, public_ref_matrix, metadata=runtime_config),
                                    config=runtime_config,
                                    optimizer_seed=seeds["optimizer_seed"],
                                ),
                            },
                            {
                                "method": "shuffled_template_n100",
                                "reference_population": "eur",
                                "template_design": "column_shuffled",
                                "reference_ids": split["public_ref"],
                                "run": lambda problem=problem, shuffled_matrix=shuffled_matrix, runtime_config=runtime_config, seeds=seeds: _with_config(
                                    reference_haplotype_reconstruction(problem, shuffled_matrix, metadata=runtime_config),
                                    config=runtime_config,
                                    optimizer_seed=seeds["optimizer_seed"],
                                ),
                            },
                        ]
                    )
                if args.run_saleem:
                    run_specs.append(
                        {
                            "method": SALEEM_RECONSTRUCTION_ATTACK,
                            "reference_population": "eur",
                            "template_design": "saleem_ld",
                            "reference_ids": split["public_ref"],
                            "run": lambda problem=problem, public_ref_matrix=public_ref_matrix, beacon_matrix=beacon_matrix, runtime_config=runtime_config, seeds=seeds: _with_config(
                                run_saleem_attack(problem, public_ref_matrix, beacon_matrix, args),
                                config=runtime_config,
                                optimizer_seed=seeds["optimizer_seed"],
                            ),
                        }
                    )
                for spec in run_specs:
                    result, wall_time_sec = run_timed(spec["run"])
                    result = attach_reference_panel_mode(result, shared=args.shared_reference_panel)
                    runtime_fields = runtime_metadata_fields(result, wall_time_sec)
                    source_overlap_count = len(set(split["beacon"]).intersection(spec["reference_ids"]))
                    overlap_rows.append(
                        {
                            "query_set_id": query_set["query_set_id"],
                            "regime": query_set["regime"],
                            "beacon_size": beacon_size,
                            "replicate_id": replicate_id,
                            "method": spec["method"],
                            "source_overlap_count": source_overlap_count,
                        }
                    )
                    evaluation = evaluate_reconstruction(
                        result,
                        truth,
                        target_counts=target_counts,
                        reference_pairwise=reference_pairwise,
                        source_sets={"beacon": split["beacon"], "template_reference": spec["reference_ids"]},
                        accepted_manifest_ids=None,
                    )
                    rows.append(
                        _record(
                            query_set=query_set,
                            beacon_size=beacon_size,
                            replicate_id=replicate_id,
                            method=spec["method"],
                            reference_population=spec["reference_population"],
                            template_design=spec["template_design"],
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
        overlap_filename=OVERLAP_FILENAME,
    )
    aggregate = _aggregate(rows)
    pairwise_delta_rows = _pairwise_deltas(rows)
    dense_delta_rows = _dense_deltas(rows)
    _write_csv(out_dir / "ld_span_query_metrics.csv", rows)
    _write_csv(out_dir / "ld_span_query_aggregate.csv", aggregate)
    _write_csv(out_dir / "ld_span_query_pairwise_deltas.csv", pairwise_delta_rows)
    if dense_delta_rows:
        _write_csv(out_dir / "ld_span_query_dense_deltas.csv", dense_delta_rows)
    _write_csv(out_dir / OVERLAP_FILENAME, overlap_rows)
    _write_json(out_dir / "ld_span_query_manifest.json", {"manifest_id": MANIFEST_ID, "config": config})
    _write_report(
        out_dir / "ld_span_query_report.md",
        aggregate_rows=aggregate,
        pairwise_rows=pairwise_delta_rows,
        dense_rows=dense_delta_rows,
        config=config,
    )
    print(f"wrote {out_dir / 'ld_span_query_report.md'}")


if __name__ == "__main__":
    main()
