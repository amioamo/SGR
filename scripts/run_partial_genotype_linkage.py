"""Partial-genotype linkage accuracy on reconstructed Beacon rows."""

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
from beacon_reconstruction.core import ExperimentConfig, BeaconAttackProblem, ReconstructionResult
from beacon_reconstruction.evaluation import align_rows_hungarian, compute_partial_genotype_linkage
from beacon_reconstruction.ld_af_baseline import SALEEM_RECONSTRUCTION_ATTACK
from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers, matrix_for_samples as _matrix_for_samples
from scripts.lib.sample_splits import load_sample_manifest, load_seed_manifest, split_samples
from scripts.lib.reference_helpers import sample_reference_ids as _sample_reference_ids, shuffle_reference_columns as _shuffle_template_columns
from scripts.ld_af_helpers import (
    add_saleem_arguments,
    run_saleem_attack,
    saleem_config_dict,
)


MANIFEST_ID = "1000g_eur_hg38_partial_genotype_linkage_v1"
MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
OUT_DIR = ROOT / "results" / "block_inheritance_v1" / "partial_genotype_linkage_shared_panel"
KNOWN_SNP_COUNTS = (10, 25, 50, 100)
LINKAGE_METHODS = (
    "matched_template_n100",
    SALEEM_RECONSTRUCTION_ATTACK,
    "pairwise_refinement",
    "independent_count",
)


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


def _run_summary_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        int(row["beacon_size"]),
        row["window_id"],
        row["replicate_id"],
        row["method"],
        int(row["known_snp_count"]),
    )


def _merge_existing_rows(
    new_rows: list[dict[str, Any]],
    *,
    metrics_path: Path,
) -> list[dict[str, Any]]:
    if not metrics_path.exists():
        return new_rows
    existing_rows = _read_csv(metrics_path)
    new_keys = {_run_summary_key(row) for row in new_rows}
    kept = [row for row in existing_rows if _run_summary_key(row) not in new_keys]
    return kept + new_rows


def _linkage_rows_for_run(
    *,
    beacon_size: int,
    window_id: str,
    replicate_id: str,
    method: str,
    beacon_ids: list[str],
    true_matrix: np.ndarray,
    reconstructed_matrix: np.ndarray,
    linkage_seed: int,
) -> list[dict[str, Any]]:
    alignment = align_rows_hungarian(true_matrix, reconstructed_matrix)
    out: list[dict[str, Any]] = []
    for known_snp_count in KNOWN_SNP_COUNTS:
        metrics = compute_partial_genotype_linkage(
            true_matrix,
            reconstructed_matrix,
            known_snp_count,
            base_seed=linkage_seed,
            alignment=alignment,
        )
        out.append(
            {
                "beacon_size": beacon_size,
                "window_id": window_id,
                "replicate_id": replicate_id,
                "method": method,
                "known_snp_count": known_snp_count,
                "n_participants": metrics.n_participants,
                "top1_accuracy": metrics.top1_accuracy,
                "top5_accuracy": metrics.top5_accuracy,
                "mean_rank_of_true_match": metrics.mean_rank_of_true_match,
                "median_rank_of_true_match": metrics.median_rank_of_true_match,
                "linkage_seed": linkage_seed,
            }
        )
    return out


def _aggregate_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    keys = sorted(
        {
            (int(row["beacon_size"]), row["method"], int(row["known_snp_count"]))
            for row in rows
        }
    )
    for beacon_size, method, known_snp_count in keys:
        selected = [
            row
            for row in rows
            if int(row["beacon_size"]) == beacon_size
            and row["method"] == method
            and int(row["known_snp_count"]) == known_snp_count
        ]
        out.append(
            {
                "beacon_size": beacon_size,
                "method": method,
                "known_snp_count": known_snp_count,
                "n_runs": len({(row["window_id"], row["replicate_id"]) for row in selected}),
                "mean_top1_accuracy": float(np.mean([float(row["top1_accuracy"]) for row in selected])),
                "mean_top5_accuracy": float(np.mean([float(row["top5_accuracy"]) for row in selected])),
                "mean_rank_of_true_match": float(
                    np.mean([float(row["mean_rank_of_true_match"]) for row in selected])
                ),
            }
        )
    return out


def _write_report(path: Path, aggregate_rows: list[dict[str, Any]], config: dict[str, Any]) -> None:
    lines = [
        "# Partial-Genotype Linkage on Reconstructed Rows",
        "",
        "## Configuration",
        "",
        f"- windows: {config['max_windows']}",
        f"- replicates: {config['max_replicates']}",
        f"- beacon sizes: {','.join(str(size) for size in config['beacon_sizes'])}",
        f"- known SNP counts: {','.join(str(size) for size in config['known_snp_counts'])}",
        f"- methods: {','.join(config['methods'])}",
        "",
        "## Aggregate Linkage Accuracy",
        "",
        "| N | method | known SNPs | runs | top-1 | top-5 | mean rank |",
        "|---:|---|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['beacon_size']} | {row['method']} | {row['known_snp_count']} | {row['n_runs']} | "
            f"{row['mean_top1_accuracy']:.4f} | {row['mean_top5_accuracy']:.4f} | "
            f"{row['mean_rank_of_true_match']:.2f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Each participant is paired to a reconstructed row after anonymous Hungarian alignment. "
            "An attacker with carrier status at k randomly selected SNPs ranks reconstructed rows by "
            "partial Hamming distance and is scored on whether the aligned row is recovered in the top ranks.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument(
        "--beacon-sizes",
        default="25,50,100,250",
        help="Comma-separated beacon sizes.",
    )
    parser.add_argument(
        "--methods",
        default=",".join(LINKAGE_METHODS),
        help="Comma-separated method ids to evaluate.",
    )
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument(
        "--merge-existing",
        action="store_true",
        help="Keep rows in out-dir CSVs that are not replaced by this run.",
    )
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--proposals-per-epoch", type=int, default=5000)
    parser.add_argument("--early-stopping", type=int, default=5)
    add_saleem_arguments(parser)
    args = parser.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    windows = _read_csv(MULTIWINDOW_DIR / "multiwindow_manifest.csv")[: max(0, args.max_windows)]
    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]
    beacon_sizes = sorted({int(value) for value in args.beacon_sizes.split(",") if value.strip()})
    methods = [value.strip() for value in args.methods.split(",") if value.strip()]
    config = {
        "max_windows": args.max_windows,
        "max_replicates": args.max_replicates,
        "beacon_sizes": beacon_sizes,
        "known_snp_counts": list(KNOWN_SNP_COUNTS),
        "methods": methods,
        "epochs": args.epochs,
        "proposals_per_epoch": args.proposals_per_epoch,
        "early_stopping": args.early_stopping,
        "template_selection_rule": TEMPLATE_SELECTION_RULE,
        **saleem_config_dict(args),
    }

    run_rows: list[dict[str, Any]] = []
    for beacon_size in beacon_sizes:
        for window in windows:
            window_id = window["window_id"]
            sample_ids, variant_ids, positions, target_alleles, carrier_matrix = _load_raw_as_target_carriers(
                Path(window["raw_path"])
            )

            for replicate_id in replicate_ids:
                seeds = seeds_by_replicate[replicate_id]
                split = split_samples(
                    sample_rows,
                    replicate_id,
                    beacon_size,
                    shared_reference_panel=True,
                )
                template_ref_ids = split["public_ref"]
                beacon_matrix = _matrix_for_samples(sample_ids, carrier_matrix, split["beacon"])
                public_ref_matrix = _matrix_for_samples(sample_ids, carrier_matrix, split["public_ref"])
                arg_ref_matrix = _matrix_for_samples(sample_ids, carrier_matrix, template_ref_ids)
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
                        metadata={"window_id": window_id, "experiment": "partial_genotype_linkage"},
                    ),
                )
                runtime_config = {
                    "window_id": window_id,
                    "replicate_id": replicate_id,
                    "beacon_size": beacon_size,
                    "config": config,
                }
                linkage_seed = int(seeds["group_selection_seed"]) + 29_003
                panel_ids = _sample_reference_ids(
                    template_ref_ids,
                    panel_size=100,
                    seed=seeds["group_selection_seed"] + 100,
                )
                panel_matrix = _matrix_for_samples(sample_ids, carrier_matrix, panel_ids)
                shuffled_matrix = _shuffle_template_columns(
                    arg_ref_matrix,
                    seed=seeds["group_selection_seed"] + 19_001,
                )

                method_results: dict[str, ReconstructionResult] = {}
                if "independent_count" in methods:
                    method_results["independent_count"] = _with_config(
                        independent_count(problem, reference_pairwise=reference_pairwise),
                        config=runtime_config,
                        optimizer_seed=seeds["optimizer_seed"],
                    )
                if "matched_template_n100" in methods:
                    method_results["matched_template_n100"] = _with_config(
                        reference_haplotype_reconstruction(problem, panel_matrix, metadata=runtime_config),
                        config=runtime_config,
                        optimizer_seed=seeds["optimizer_seed"],
                    )
                if "pairwise_refinement" in methods:
                    method_results["pairwise_refinement"] = _with_config(
                        pairwise_refinement(
                            problem,
                            reference_pairwise,
                            epochs=args.epochs,
                            proposals_per_epoch=args.proposals_per_epoch,
                            early_stopping=args.early_stopping,
                        ),
                        config=runtime_config,
                        optimizer_seed=seeds["optimizer_seed"],
                    )
                if "shuffled_template_n100" in methods:
                    method_results["shuffled_template_n100"] = _with_config(
                        reference_haplotype_reconstruction(problem, shuffled_matrix, metadata=runtime_config),
                        config=runtime_config,
                        optimizer_seed=seeds["optimizer_seed"],
                    )
                if SALEEM_RECONSTRUCTION_ATTACK in methods and args.run_saleem:
                    method_results[SALEEM_RECONSTRUCTION_ATTACK] = _with_config(
                        run_saleem_attack(problem, public_ref_matrix, beacon_matrix, args),
                        config=runtime_config,
                        optimizer_seed=seeds["optimizer_seed"],
                    )

                for method in methods:
                    if method not in method_results:
                        continue
                    run_rows.extend(
                        _linkage_rows_for_run(
                            beacon_size=beacon_size,
                            window_id=window_id,
                            replicate_id=replicate_id,
                            method=method,
                            beacon_ids=split["beacon"],
                            true_matrix=beacon_matrix,
                            reconstructed_matrix=method_results[method].reconstructed_binary,
                            linkage_seed=linkage_seed,
                        )
                    )
                print(f"completed {split_id}")

    if args.merge_existing:
        run_rows = _merge_existing_rows(
            run_rows,
            metrics_path=out_dir / "partial_genotype_linkage_metrics.csv",
        )
    aggregate = _aggregate_rows(run_rows)
    _write_csv(out_dir / "partial_genotype_linkage_metrics.csv", run_rows)
    _write_csv(out_dir / "partial_genotype_linkage_aggregate.csv", aggregate)
    _write_json(out_dir / "partial_genotype_linkage_manifest.json", {"manifest_id": MANIFEST_ID, "config": config})
    _write_report(out_dir / "partial_genotype_linkage_report.md", aggregate, config)
    print(f"wrote {out_dir / 'partial_genotype_linkage_report.md'}")


if __name__ == "__main__":
    main()
