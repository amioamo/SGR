"""Run cross-ancestry reference matching for genomic data-sharing privacy risk."""

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
from beacon_reconstruction.evaluation import EvaluationResult, evaluate_reconstruction
from scripts.ld_af_helpers import (
    add_saleem_arguments,
    merge_saleem_rows,
    run_saleem_attack,
    saleem_config_dict,
)
from scripts.lib.multiwindow_data import load_raw_carriers
from scripts.experiment_runtime_helpers import (
    attach_reference_panel_mode,
    run_timed,
    runtime_metadata_fields,
)


MANIFEST_ID = "1000g_hg38_multiancestry_reference_matching_v1"
DEFAULT_PANEL_DIR = ROOT / "data" / "multiwindow" / "ancestry_panels"
DEFAULT_OUT_DIR = ROOT / "results" / "experiments" / "multiancestry_reference_matching"
DEFAULT_SHARED_OUT_DIR = (
    ROOT / "results" / "experiments" / "multiancestry_reference_matching_shared_panel"
)
DEFAULT_MATCHED_FRACTIONS = (0, 25, 50, 75, 100)
DEFAULT_BEACON_SIZES = (25, 50, 100, 250)
SALEEM_MIXED_PREFIX = "saleem_mixed_matched_"
LEGACY_SALEEM_METHODS = ("saleem_matched_reconstruction", "saleem_wrong_reconstruction")


def saleem_methods_for_fractions(matched_fractions: list[int]) -> tuple[str, ...]:
    mixed = tuple(f"{SALEEM_MIXED_PREFIX}{fraction}" for fraction in matched_fractions)
    return LEGACY_SALEEM_METHODS + mixed


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


def _matrix_for_indices(matrix: np.ndarray, indices: np.ndarray) -> np.ndarray:
    return matrix[indices.astype(int), :]


def _select_indices(total: int, count: int, *, seed: int, exclude: set[int] | None = None) -> np.ndarray:
    if exclude is None:
        exclude = set()
    candidates = np.asarray([index for index in range(total) if index not in exclude], dtype=int)
    if len(candidates) < count:
        raise ValueError(f"Need {count} samples but only {len(candidates)} are available.")
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(candidates, size=count, replace=False).astype(int))


def _split_target_indices(
    total: int,
    *,
    beacon_size: int,
    seed: int,
    shared_reference_panel: bool = False,
) -> dict[str, np.ndarray]:
    ref_count = 100
    required = beacon_size + ref_count if shared_reference_panel else beacon_size + 2 * ref_count
    if total < required:
        raise ValueError(f"Need {required} target samples for N={beacon_size}, found {total}.")
    rng = np.random.default_rng(seed)
    selected = rng.permutation(total)[:required].astype(int)
    public_ref = np.sort(selected[beacon_size : beacon_size + ref_count])
    if shared_reference_panel:
        template_ref = public_ref.copy()
    else:
        template_ref = np.sort(selected[beacon_size + ref_count : required])
    return {
        "beacon": np.sort(selected[:beacon_size]),
        "public_ref": public_ref,
        "template_ref": template_ref,
    }


def _sample_template_panel(matrix: np.ndarray, indices: np.ndarray, count: int, *, seed: int) -> tuple[np.ndarray, np.ndarray]:
    selected = _select_indices(len(indices), count, seed=seed)
    panel_indices = indices[selected]
    return panel_indices, _matrix_for_indices(matrix, panel_indices)


def _mixed_template_panel(
    *,
    matched_matrix: np.ndarray,
    matched_indices: np.ndarray,
    wrong_matrix: np.ndarray,
    wrong_indices: np.ndarray,
    matched_fraction: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    matched_count = int(round(100 * matched_fraction / 100.0))
    wrong_count = 100 - matched_count
    matched_selected = np.array([], dtype=int)
    wrong_selected = np.array([], dtype=int)
    parts: list[np.ndarray] = []
    if matched_count:
        matched_selected, matched_part = _sample_template_panel(matched_matrix, matched_indices, matched_count, seed=seed + 11)
        parts.append(matched_part)
    if wrong_count:
        wrong_selected, wrong_part = _sample_template_panel(wrong_matrix, wrong_indices, wrong_count, seed=seed + 29)
        parts.append(wrong_part)
    return np.concatenate([matched_selected, wrong_selected]), np.vstack(parts)


def _with_config(result: ReconstructionResult, *, config: dict[str, Any], optimizer_seed: int) -> ReconstructionResult:
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
    target_ancestry: str,
    reference_ancestry: str,
    beacon_size: int,
    window_id: str,
    replicate_id: int,
    method: str,
    reference_design: str,
    matched_fraction: int | None,
    result: ReconstructionResult,
    evaluation: EvaluationResult,
    source_overlap_count: int,
    runtime_fields: dict[str, Any] | None = None,
) -> dict[str, Any]:
    metrics = evaluation.metrics
    row = {
        "target_ancestry": target_ancestry,
        "reference_ancestry": reference_ancestry,
        "beacon_size": beacon_size,
        "window_id": window_id,
        "replicate_id": replicate_id,
        "method": method,
        "reference_design": reference_design,
        "matched_fraction": "" if matched_fraction is None else matched_fraction,
        "precision": metrics.precision,
        "recall": metrics.recall,
        "f1": metrics.f1,
        "pairwise_error": metrics.pairwise_error,
        "alignment_total_hamming_cost": metrics.alignment_total_hamming_cost,
        "max_count_error": metrics.max_count_error,
        "source_overlap_count": source_overlap_count,
        "invalid_run": evaluation.invalid_run.is_invalid,
        "invalid_reasons": ";".join(evaluation.invalid_run.reasons),
    }
    if runtime_fields:
        row.update(runtime_fields)
    return row


def _aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    keys = sorted(
        {
            (
                row["target_ancestry"],
                row["reference_ancestry"],
                int(row["beacon_size"]),
                row["method"],
                row["reference_design"],
                str(row["matched_fraction"]),
            )
            for row in rows
        }
    )
    for target, reference, beacon_size, method, design, fraction in keys:
        selected = [
            row
            for row in rows
            if row["target_ancestry"] == target
            and row["reference_ancestry"] == reference
            and int(row["beacon_size"]) == beacon_size
            and row["method"] == method
            and row["reference_design"] == design
            and str(row["matched_fraction"]) == fraction
        ]
        out.append(
            {
                "target_ancestry": target,
                "reference_ancestry": reference,
                "beacon_size": beacon_size,
                "method": method,
                "reference_design": design,
                "matched_fraction": fraction,
                "n_runs": len(selected),
                "mean_f1": float(np.mean([float(row["f1"]) for row in selected])),
                "median_f1": float(np.median([float(row["f1"]) for row in selected])),
                "mean_hamming": float(np.mean([float(row["alignment_total_hamming_cost"]) for row in selected])),
                "source_overlap_count": int(sum(int(row["source_overlap_count"]) for row in selected)),
                "invalid_run_count": int(sum(bool(row["invalid_run"]) for row in selected)),
            }
        )
        if "wall_time_sec" in selected[0]:
            wall_times = [float(row["wall_time_sec"]) for row in selected]
            out[-1]["mean_wall_time_sec"] = float(np.mean(wall_times))
            out[-1]["median_wall_time_sec"] = float(np.median(wall_times))
    return out


def _write_report(path: Path, *, aggregate_rows: list[dict[str, Any]], config: dict[str, Any]) -> None:
    lines = [
        "# Multi-Ancestry Reference Matching",
        "",
        "## Configuration",
        "",
        f"- target ancestries: {','.join(config['target_ancestries'])}",
        f"- beacon sizes: {','.join(str(size) for size in config['beacon_sizes'])}",
        f"- windows: {config['max_windows']}",
        f"- replicates: {config['replicates']}",
        f"- matched fractions: {','.join(str(value) for value in config['matched_fractions'])}",
        f"- reference panel mode: {config.get('reference_panel_mode', 'dual_disjoint')}",
        f"- pairwise budget: {config['epochs']} epochs x {config['proposals_per_epoch']} proposals",
        "",
        "## Aggregate Metrics",
        "",
        "| target | reference | N | method | design | matched_fraction | n_runs | mean_f1 | median_f1 | overlap | invalid |",
        "|---|---|---:|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['target_ancestry']} | {row['reference_ancestry']} | {row['beacon_size']} | "
            f"{row['method']} | {row['reference_design']} | {row['matched_fraction']} | {row['n_runs']} | "
            f"{row['mean_f1']:.6f} | {row['median_f1']:.6f} | {row['source_overlap_count']} | "
            f"{row['invalid_run_count']} |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-dir", type=Path, default=DEFAULT_PANEL_DIR)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--target-ancestries", default="EUR,AFR,EAS,SAS,AMR")
    parser.add_argument("--beacon-sizes", default="25,50,100,250")
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--matched-fractions", default="0,25,50,75,100")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--proposals-per-epoch", type=int, default=5000)
    parser.add_argument("--early-stopping", type=int, default=5)
    parser.add_argument(
        "--shared-reference-panel",
        action="store_true",
        help="Use one R=100 reference panel for LD/Saleem/pairwise and template rows (disjoint from beacon only).",
    )
    add_saleem_arguments(parser)
    args = parser.parse_args()

    out_dir = args.out_dir or (DEFAULT_SHARED_OUT_DIR if args.shared_reference_panel else DEFAULT_OUT_DIR)
    args.out_dir = out_dir

    panel_manifest = _read_csv(args.panel_dir / "ancestry_panel_manifest.csv")
    by_key = {(row["ancestry"], row["window_id"]): Path(row["raw_path"]) for row in panel_manifest}
    windows = sorted({row["window_id"] for row in panel_manifest})[: max(0, args.max_windows)]
    target_ancestries = [value.strip().upper() for value in args.target_ancestries.split(",") if value.strip()]
    beacon_sizes = sorted({int(value) for value in args.beacon_sizes.split(",") if value.strip()})
    matched_fractions = sorted({int(value) for value in args.matched_fractions.split(",") if value.strip()})
    saleem_methods = saleem_methods_for_fractions(matched_fractions)
    config = {
        "target_ancestries": target_ancestries,
        "beacon_sizes": beacon_sizes,
        "max_windows": args.max_windows,
        "replicates": args.replicates,
        "matched_fractions": matched_fractions,
        "epochs": args.epochs,
        "proposals_per_epoch": args.proposals_per_epoch,
        "early_stopping": args.early_stopping,
        "panel_dir": str(args.panel_dir),
        "out_dir": str(args.out_dir),
        **saleem_config_dict(args),
        "saleem_mixed_fractions": matched_fractions,
        "reference_panel_mode": "shared" if args.shared_reference_panel else "dual_disjoint",
    }

    wrong_by_target = {ancestry: ("AFR" if ancestry == "EUR" else "EUR") for ancestry in target_ancestries}
    rows: list[dict[str, Any]] = []
    for window_id in windows:
        for target_ancestry in target_ancestries:
            wrong_ancestry = wrong_by_target[target_ancestry]
            target_ids, variant_ids, positions, target_alleles, target_matrix = load_raw_carriers(
                by_key[(target_ancestry, window_id)]
            )
            wrong_ids, wrong_variant_ids, _wrong_positions, _wrong_targets, wrong_matrix = load_raw_carriers(
                by_key[(wrong_ancestry, window_id)],
                target_alleles=target_alleles,
            )
            if wrong_variant_ids != variant_ids:
                raise ValueError(f"Variant mismatch for {window_id}: {target_ancestry} vs {wrong_ancestry}.")
            wrong_public_indices = _select_indices(len(wrong_ids), 100, seed=910_000 + len(window_id) + len(wrong_ancestry))
            if args.shared_reference_panel:
                wrong_template_indices = wrong_public_indices.copy()
            else:
                wrong_template_indices = _select_indices(
                    len(wrong_ids),
                    100,
                    seed=920_000 + len(window_id) + len(wrong_ancestry),
                    exclude=set(int(index) for index in wrong_public_indices),
                )
            wrong_public_matrix = _matrix_for_indices(wrong_matrix, wrong_public_indices)
            wrong_template_matrix = _matrix_for_indices(wrong_matrix, wrong_template_indices)

            for beacon_size in beacon_sizes:
                for replicate_id in range(args.replicates):
                    seed_base = 100_000 + 10_000 * replicate_id + 100 * beacon_size + len(window_id)
                    split = _split_target_indices(
                        len(target_ids),
                        beacon_size=beacon_size,
                        seed=seed_base,
                        shared_reference_panel=args.shared_reference_panel,
                    )
                    beacon_matrix = _matrix_for_indices(target_matrix, split["beacon"])
                    matched_public_matrix = _matrix_for_indices(target_matrix, split["public_ref"])
                    matched_template_matrix = _matrix_for_indices(target_matrix, split["template_ref"])
                    target_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                    matched_pairwise = pairwise_moments(matched_public_matrix)
                    wrong_pairwise = pairwise_moments(wrong_public_matrix)
                    split_id = f"{window_id}_{target_ancestry}_replicate_{replicate_id}_N{beacon_size}"
                    runtime_config = {
                        "window_id": window_id,
                        "target_ancestry": target_ancestry,
                        "wrong_ancestry": wrong_ancestry,
                        "beacon_size": beacon_size,
                        "replicate_id": replicate_id,
                    }
                    problem = BeaconAttackProblem(
                        manifest_id=MANIFEST_ID,
                        split_id=split_id,
                        beacon_size=beacon_size,
                        variant_ids=variant_ids,
                        positions=positions,
                        query_responses=target_counts,
                        target_counts=target_counts,
                        public_ld=matched_pairwise,
                        config=ExperimentConfig(
                            initialization_seed=seed_base + 1,
                            optimizer_seed=seed_base + 2,
                            group_selection_seed=seed_base + 3,
                            metadata={"window_id": window_id, "experiment": "multiancestry_reference_matching"},
                        ),
                    )
                    beacon_sample_ids = [target_ids[int(index)] for index in split["beacon"]]
                    truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(beacon_sample_ids))

                    run_spec_defs: list[dict[str, Any]] = []
                    if not args.saleem_only:
                        run_spec_defs.extend(
                            [
                                {
                                    "method": "independent_count",
                                    "reference_ancestry": "none",
                                    "reference_design": "count_only",
                                    "matched_fraction": None,
                                    "reference_sample_ids": [],
                                    "run": lambda problem=problem, matched_pairwise=matched_pairwise, runtime_config=runtime_config, seed_base=seed_base: _with_config(
                                        independent_count(problem, reference_pairwise=matched_pairwise),
                                        config=runtime_config,
                                        optimizer_seed=seed_base + 2,
                                    ),
                                },
                                {
                                    "method": "pairwise_matched_refinement",
                                    "reference_ancestry": target_ancestry,
                                    "reference_design": "matched_pairwise",
                                    "matched_fraction": 100,
                                    "reference_sample_ids": [
                                        target_ids[int(index)] for index in split["public_ref"]
                                    ],
                                    "run": lambda problem=problem, matched_pairwise=matched_pairwise, runtime_config=runtime_config, seed_base=seed_base: _with_config(
                                        pairwise_refinement(
                                            problem,
                                            matched_pairwise,
                                            epochs=args.epochs,
                                            proposals_per_epoch=args.proposals_per_epoch,
                                            early_stopping=args.early_stopping,
                                        ),
                                        config=runtime_config,
                                        optimizer_seed=seed_base + 2,
                                    ),
                                },
                                {
                                    "method": "pairwise_wrong_refinement",
                                    "reference_ancestry": wrong_ancestry,
                                    "reference_design": "wrong_pairwise",
                                    "matched_fraction": 0,
                                    "reference_sample_ids": [
                                        wrong_ids[int(index)] for index in wrong_public_indices
                                    ],
                                    "run": lambda problem=problem, wrong_pairwise=wrong_pairwise, runtime_config=runtime_config, seed_base=seed_base: _with_config(
                                        pairwise_refinement(
                                            problem,
                                            wrong_pairwise,
                                            epochs=args.epochs,
                                            proposals_per_epoch=args.proposals_per_epoch,
                                            early_stopping=args.early_stopping,
                                        ),
                                        config=runtime_config,
                                        optimizer_seed=seed_base + 2,
                                    ),
                                },
                            ]
                        )
                    if args.run_saleem:
                        for fraction in matched_fractions:
                            _selected_indices, mixed_public_matrix = _mixed_template_panel(
                                matched_matrix=target_matrix,
                                matched_indices=split["public_ref"],
                                wrong_matrix=wrong_matrix,
                                wrong_indices=wrong_public_indices,
                                matched_fraction=fraction,
                                seed=seed_base + 3 + fraction + 1000,
                            )
                            mixed_pairwise = pairwise_moments(mixed_public_matrix)
                            mixed_reference_ids: list[str] = []
                            if fraction:
                                mixed_reference_ids.extend(
                                    [target_ids[int(index)] for index in split["public_ref"]]
                                )
                            if fraction < 100:
                                mixed_reference_ids.extend(
                                    [wrong_ids[int(index)] for index in wrong_public_indices]
                                )
                            mixed_problem = BeaconAttackProblem(
                                manifest_id=MANIFEST_ID,
                                split_id=split_id,
                                beacon_size=beacon_size,
                                variant_ids=variant_ids,
                                positions=positions,
                                query_responses=target_counts,
                                target_counts=target_counts,
                                public_ld=mixed_pairwise,
                                config=problem.config,
                            )
                            run_spec_defs.append(
                                {
                                    "method": f"{SALEEM_MIXED_PREFIX}{fraction}",
                                    "reference_ancestry": f"{target_ancestry}+{wrong_ancestry}",
                                    "reference_design": "mixed_saleem",
                                    "matched_fraction": fraction,
                                    "reference_sample_ids": mixed_reference_ids,
                                    "run": lambda mixed_problem=mixed_problem, mixed_public_matrix=mixed_public_matrix, beacon_matrix=beacon_matrix, runtime_config=runtime_config, seed_base=seed_base: _with_config(
                                        run_saleem_attack(
                                            mixed_problem,
                                            mixed_public_matrix,
                                            beacon_matrix,
                                            args,
                                        ),
                                        config=runtime_config,
                                        optimizer_seed=seed_base + 2,
                                    ),
                                }
                            )

                    if not args.saleem_only:
                        for fraction in matched_fractions:
                            _selected_indices, mixed_matrix = _mixed_template_panel(
                                matched_matrix=target_matrix,
                                matched_indices=split["template_ref"],
                                wrong_matrix=wrong_matrix,
                                wrong_indices=wrong_template_indices,
                                matched_fraction=fraction,
                                seed=seed_base + 3 + fraction,
                            )
                            reference_ids = []
                            if fraction:
                                reference_ids.extend(
                                    [target_ids[int(index)] for index in split["template_ref"]]
                                )
                            if fraction < 100:
                                reference_ids.extend(
                                    [wrong_ids[int(index)] for index in wrong_template_indices]
                                )
                            run_spec_defs.append(
                                {
                                    "method": f"template_mixed_matched_{fraction}",
                                    "reference_ancestry": f"{target_ancestry}+{wrong_ancestry}",
                                    "reference_design": "mixed_template",
                                    "matched_fraction": fraction,
                                    "reference_sample_ids": reference_ids,
                                    "run": lambda problem=problem, mixed_matrix=mixed_matrix, runtime_config=runtime_config, seed_base=seed_base: _with_config(
                                        reference_haplotype_reconstruction(problem, mixed_matrix, metadata=runtime_config),
                                        config=runtime_config,
                                        optimizer_seed=seed_base + 2,
                                    ),
                                }
                            )

                    for spec_def in run_spec_defs:
                        result, wall_time_sec = run_timed(spec_def["run"])
                        result = attach_reference_panel_mode(result, shared=args.shared_reference_panel)
                        runtime_fields = runtime_metadata_fields(result, wall_time_sec)
                        source_overlap_count = len(
                            set(beacon_sample_ids).intersection(spec_def["reference_sample_ids"])
                        )
                        evaluation = evaluate_reconstruction(
                            result,
                            truth,
                            target_counts=target_counts,
                            reference_pairwise=matched_pairwise,
                            source_sets={
                                "beacon": beacon_sample_ids,
                                "reference": spec_def["reference_sample_ids"],
                            },
                            accepted_manifest_ids=None,
                        )
                        rows.append(
                            _record(
                                target_ancestry=target_ancestry,
                                reference_ancestry=spec_def["reference_ancestry"],
                                beacon_size=beacon_size,
                                window_id=window_id,
                                replicate_id=replicate_id,
                                method=spec_def["method"],
                                reference_design=spec_def["reference_design"],
                                matched_fraction=spec_def["matched_fraction"],
                                result=result,
                                evaluation=evaluation,
                                source_overlap_count=source_overlap_count,
                                runtime_fields=runtime_fields,
                            )
                        )
                    print(f"completed {split_id}")

    rows = merge_saleem_rows(
        new_rows=rows,
        merge_metrics=args.merge_metrics,
        saleem_only=args.saleem_only,
        saleem_methods=saleem_methods,
    )
    aggregate = _aggregate(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(args.out_dir / "multiancestry_reference_matching_metrics.csv", rows)
    _write_csv(args.out_dir / "multiancestry_reference_matching_aggregate.csv", aggregate)
    _write_json(args.out_dir / "multiancestry_reference_matching_manifest.json", {"manifest_id": MANIFEST_ID, "config": config})
    _write_report(args.out_dir / "multiancestry_reference_matching_report.md", aggregate_rows=aggregate, config=config)
    print(f"wrote {args.out_dir / 'multiancestry_reference_matching_report.md'}")


if __name__ == "__main__":
    main()
