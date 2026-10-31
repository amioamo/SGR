"""Evaluate count-release perturbations for genomic data-sharing privacy risk."""

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
from scripts.experiment_runtime_helpers import (
    attach_reference_panel_mode,
    run_timed,
    runtime_metadata_fields,
)


MANIFEST_ID = "1000g_eur_hg38_count_release_perturbation_v1"
DEFAULT_MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
DEFAULT_OUT_DIR = ROOT / "results" / "experiments" / "defended_beacon_extended"
DEFAULT_SHARED_OUT_DIR = ROOT / "results" / "experiments" / "defended_beacon_extended_shared_panel"
DEFAULT_PROFILES = (
    "exact_count",
    "round2",
    "suppress_le2",
    "round2_suppress_le2",
    "noise_pm1_suppress_le2",
    "frequency_round_0.05",
    "query_limit_100",
    "query_limit_250",
    "min_cohort_50",
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


def _all_indices(n_variants: int) -> np.ndarray:
    return np.arange(n_variants, dtype=int)


def _query_limit_indices(n_variants: int, limit: int, *, seed: int) -> np.ndarray:
    if limit >= n_variants:
        return _all_indices(n_variants)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_variants, size=max(0, limit), replace=False).astype(int))


def _profile_indices(
    profile: str,
    n_variants: int,
    *,
    beacon_size: int,
    seed: int,
) -> np.ndarray:
    if profile.startswith("query_limit_"):
        return _query_limit_indices(n_variants, int(profile.rsplit("_", 1)[1]), seed=seed)
    if profile == "min_cohort_50" and beacon_size < 50:
        return np.array([], dtype=int)
    return _all_indices(n_variants)


def _apply_profile(counts: np.ndarray, *, profile: str, beacon_size: int, seed: int) -> np.ndarray:
    released = np.asarray(counts, dtype=np.int64).copy()
    if profile in {"exact_count", "min_cohort_50"} or profile.startswith("query_limit_"):
        pass
    elif profile == "round2":
        released = ((released + 1) // 2) * 2
    elif profile == "suppress_le2":
        released = np.where(released <= 2, 0, released)
    elif profile == "round2_suppress_le2":
        released = np.where(released <= 2, 0, released)
        released = ((released + 1) // 2) * 2
    elif profile == "noise_pm1_suppress_le2":
        rng = np.random.default_rng(seed)
        noise = rng.choice(np.asarray([-1, 0, 1], dtype=np.int64), size=len(released), p=[0.25, 0.5, 0.25])
        released = np.where(released <= 2, 0, released + noise)
    elif profile.startswith("frequency_round_"):
        step = float(profile.rsplit("_", 1)[1])
        released_freq = np.round((released / float(beacon_size)) / step) * step
        released = np.rint(released_freq * float(beacon_size)).astype(np.int64)
    else:
        raise ValueError(f"Unsupported count-release profile: {profile}.")
    return np.clip(released, 0, beacon_size).astype(np.int64, copy=False)


def _with_metadata(result: ReconstructionResult, *, config: dict[str, Any], optimizer_seed: int) -> ReconstructionResult:
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


def _release_stats(released_counts: np.ndarray, true_counts: np.ndarray, *, released_variant_count: int, original_variant_count: int) -> dict[str, float | int]:
    delta = np.asarray(released_counts, dtype=np.int64) - np.asarray(true_counts, dtype=np.int64)
    return {
        "released_variant_count": int(released_variant_count),
        "original_variant_count": int(original_variant_count),
        "released_variant_fraction": float(released_variant_count / original_variant_count) if original_variant_count else float("nan"),
        "count_mae_to_true": float(np.mean(np.abs(delta))) if delta.size else float("nan"),
        "count_max_abs_error": int(np.max(np.abs(delta), initial=0)),
        "changed_count_fraction": float(np.mean(delta != 0)) if delta.size else float("nan"),
    }


def _record(
    *,
    beacon_size: int,
    window_id: str,
    replicate_id: str,
    profile: str,
    method: str,
    result: ReconstructionResult,
    evaluation: EvaluationResult,
    release_stats: dict[str, float | int],
    source_overlap_count: int,
    runtime_fields: dict[str, Any] | None = None,
) -> dict[str, Any]:
    metrics = evaluation.metrics
    row = {
        "beacon_size": beacon_size,
        "window_id": window_id,
        "replicate_id": replicate_id,
        "release_profile": profile,
        "method": method,
        "method_name": result.method_name,
        "precision": metrics.precision,
        "recall": metrics.recall,
        "f1": metrics.f1,
        "pairwise_error": metrics.pairwise_error,
        "alignment_total_hamming_cost": metrics.alignment_total_hamming_cost,
        "max_count_error_to_released_counts": metrics.max_count_error,
        "source_overlap_count": source_overlap_count,
        **release_stats,
        "invalid_run": evaluation.invalid_run.is_invalid,
        "invalid_reasons": ";".join(evaluation.invalid_run.reasons),
    }
    if runtime_fields:
        row.update(runtime_fields)
    return row


def _aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for beacon_size, profile, method in sorted({(int(row["beacon_size"]), row["release_profile"], row["method"]) for row in rows}):
        selected = [
            row
            for row in rows
            if int(row["beacon_size"]) == beacon_size and row["release_profile"] == profile and row["method"] == method
        ]
        out.append(
            {
                "beacon_size": beacon_size,
                "release_profile": profile,
                "method": method,
                "n_runs": len(selected),
                "mean_f1": float(np.mean([float(row["f1"]) for row in selected])),
                "median_f1": float(np.median([float(row["f1"]) for row in selected])),
                "mean_hamming": float(np.mean([float(row["alignment_total_hamming_cost"]) for row in selected])),
                "mean_released_variant_fraction": float(np.mean([float(row["released_variant_fraction"]) for row in selected])),
                "mean_count_mae_to_true": float(np.mean([float(row["count_mae_to_true"]) for row in selected])),
                "mean_changed_count_fraction": float(np.mean([float(row["changed_count_fraction"]) for row in selected])),
                "total_source_overlap_count": int(sum(int(row["source_overlap_count"]) for row in selected)),
                "invalid_run_count": int(sum(bool(row["invalid_run"]) for row in selected)),
            }
        )
        if "wall_time_sec" in selected[0]:
            wall_times = [float(row["wall_time_sec"]) for row in selected]
            out[-1]["mean_wall_time_sec"] = float(np.mean(wall_times))
    return out


def _metric_row_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        int(row["beacon_size"]),
        row["window_id"],
        row["replicate_id"],
        row["release_profile"],
        row["method"],
    )


def _merge_existing_metrics(
    new_rows: list[dict[str, Any]],
    *,
    metrics_path: Path,
) -> list[dict[str, Any]]:
    if not metrics_path.exists():
        return new_rows
    existing_rows = _read_csv(metrics_path)
    new_keys = {_metric_row_key(row) for row in new_rows}
    kept = [row for row in existing_rows if _metric_row_key(row) not in new_keys]
    return kept + new_rows


def _profile_deltas(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {
        (int(row["beacon_size"]), row["window_id"], row["replicate_id"], row["release_profile"], row["method"]): row
        for row in rows
    }
    units = sorted({(int(row["beacon_size"]), row["window_id"], row["replicate_id"]) for row in rows})
    methods = sorted({row["method"] for row in rows})
    profiles = sorted({row["release_profile"] for row in rows} - {"exact_count"})
    out: list[dict[str, Any]] = []
    for beacon_size, window_id, replicate_id in units:
        for profile in profiles:
            for method in methods:
                exact_key = (beacon_size, window_id, replicate_id, "exact_count", method)
                profile_key = (beacon_size, window_id, replicate_id, profile, method)
                if exact_key not in by_key or profile_key not in by_key:
                    continue
                exact = by_key[exact_key]
                perturbed = by_key[profile_key]
                out.append(
                    {
                        "beacon_size": beacon_size,
                        "window_id": window_id,
                        "replicate_id": replicate_id,
                        "release_profile": profile,
                        "method": method,
                        "delta_f1_vs_exact_count": float(perturbed["f1"]) - float(exact["f1"]),
                        "delta_hamming_vs_exact_count": float(perturbed["alignment_total_hamming_cost"])
                        - float(exact["alignment_total_hamming_cost"]),
                    }
                )
    return out


def _skipped_row(*, beacon_size: int, window_id: str, replicate_id: str, profile: str, reason: str, original_variant_count: int) -> dict[str, Any]:
    return {
        "beacon_size": beacon_size,
        "window_id": window_id,
        "replicate_id": replicate_id,
        "release_profile": profile,
        "reason": reason,
        "released_variant_count": 0,
        "original_variant_count": original_variant_count,
    }


def _empty_release_reason(profile: str) -> str:
    if profile == "min_cohort_50":
        return "minimum cohort size rule released no variant counts"
    return "release profile released no variant counts"


def _variant_utility_rows(
    *,
    beacon_size: int,
    window_id: str,
    replicate_id: str,
    profile: str,
    variant_ids: list[str],
    true_counts_full: np.ndarray,
    released_indices: np.ndarray,
    released_counts: np.ndarray,
) -> list[dict[str, Any]]:
    released_lookup = {int(index): int(released_counts[pos]) for pos, index in enumerate(released_indices)}
    out: list[dict[str, Any]] = []
    for index, variant_id in enumerate(variant_ids):
        released = index in released_lookup
        released_count = released_lookup.get(index)
        true_count = int(true_counts_full[index])
        abs_error = None if released_count is None else abs(released_count - true_count)
        out.append(
            {
                "beacon_size": beacon_size,
                "window_id": window_id,
                "replicate_id": replicate_id,
                "release_profile": profile,
                "variant_index": index,
                "variant_id": variant_id,
                "released": released,
                "true_count": true_count,
                "released_count": "" if released_count is None else released_count,
                "absolute_count_error": "" if abs_error is None else abs_error,
                "true_frequency": true_count / float(beacon_size),
                "released_frequency": "" if released_count is None else released_count / float(beacon_size),
            }
        )
    return out


def _write_report(path: Path, *, aggregate_rows: list[dict[str, Any]], delta_rows: list[dict[str, Any]], skipped_rows: list[dict[str, Any]], config: dict[str, Any]) -> None:
    lines = [
        "# Count-Release Perturbation Extensions",
        "",
        "## Configuration",
        "",
        f"- windows: {config['max_windows']}",
        f"- replicates: {config['max_replicates']}",
        f"- beacon sizes: {','.join(str(size) for size in config['beacon_sizes'])}",
        f"- profiles: {','.join(config['release_profiles'])}",
        f"- reference panel mode: {config.get('reference_panel_mode', 'dual_disjoint')}",
        f"- pairwise budget: {config['epochs']} epochs x {config['proposals_per_epoch']} proposals",
        "",
        "## Aggregate Metrics",
        "",
        "| N | profile | method | n_runs | mean_f1 | median_f1 | released_variant_fraction | count_mae_to_true | changed_count_fraction | invalid |",
        "|---:|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['beacon_size']} | {row['release_profile']} | {row['method']} | {row['n_runs']} | "
            f"{row['mean_f1']:.6f} | {row['median_f1']:.6f} | {row['mean_released_variant_fraction']:.6f} | "
            f"{row['mean_count_mae_to_true']:.6f} | {row['mean_changed_count_fraction']:.6f} | {row['invalid_run_count']} |"
        )

    lines.extend(
        [
            "",
            "## Privacy-Risk Change Versus Exact Counts",
            "",
            "| N | profile | method | mean_delta_f1 | median_delta_f1 | wins | ties | losses |",
            "|---:|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for beacon_size, profile, method in sorted({(int(row["beacon_size"]), row["release_profile"], row["method"]) for row in delta_rows}):
        selected = [
            row
            for row in delta_rows
            if int(row["beacon_size"]) == beacon_size and row["release_profile"] == profile and row["method"] == method
        ]
        deltas = np.asarray([float(row["delta_f1_vs_exact_count"]) for row in selected], dtype=float)
        lines.append(
            f"| {beacon_size} | {profile} | {method} | {float(np.mean(deltas)):.6f} | {float(np.median(deltas)):.6f} | "
            f"{int(np.sum(deltas > 0))} | {int(np.sum(deltas == 0))} | {int(np.sum(deltas < 0))} |"
        )

    if skipped_rows:
        lines.extend(["", "## Skipped Release Profiles", ""])
        for row in skipped_rows:
            lines.append(
                f"- N={row['beacon_size']} {row['window_id']} replicate {row['replicate_id']} "
                f"{row['release_profile']}: {row['reason']}."
            )

    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "This result separates privacy-risk reduction from utility loss. Utility loss is reported as count distortion "
            "and released-variant fraction; participant reconstruction quality is reported separately.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--beacon-sizes", default="25,50")
    parser.add_argument("--profiles", default=",".join(DEFAULT_PROFILES))
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--proposals-per-epoch", type=int, default=5000)
    parser.add_argument("--early-stopping", type=int, default=5)
    parser.add_argument(
        "--shared-reference-panel",
        action="store_true",
        help="Use one R=100 reference panel for LD/Saleem/pairwise and template rows (disjoint from beacon only).",
    )
    parser.add_argument("--multiwindow-dir", type=Path, default=DEFAULT_MULTIWINDOW_DIR)
    parser.add_argument("--out-dir", type=Path, default=None)
    add_saleem_arguments(parser)
    parser.add_argument(
        "--merge-existing",
        action="store_true",
        help="Keep rows in out-dir/count_release_perturbation_metrics.csv that are not replaced by this run.",
    )
    args = parser.parse_args()

    args.out_dir = args.out_dir or (DEFAULT_SHARED_OUT_DIR if args.shared_reference_panel else DEFAULT_OUT_DIR)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    windows = _read_csv(args.multiwindow_dir / "multiwindow_manifest.csv")[: max(0, args.max_windows)]
    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]
    beacon_sizes = sorted({int(value) for value in args.beacon_sizes.split(",") if value.strip()})
    profiles = [value.strip() for value in args.profiles.split(",") if value.strip()]
    config = {
        "max_windows": args.max_windows,
        "max_replicates": args.max_replicates,
        "beacon_sizes": beacon_sizes,
        "release_profiles": profiles,
        "epochs": args.epochs,
        "proposals_per_epoch": args.proposals_per_epoch,
        "early_stopping": args.early_stopping,
        "multiwindow_dir": str(args.multiwindow_dir),
        "out_dir": str(args.out_dir),
        **saleem_config_dict(args),
        "reference_panel_mode": "shared" if args.shared_reference_panel else "dual_disjoint",
    }

    rows: list[dict[str, Any]] = []
    utility_rows: list[dict[str, Any]] = []
    skipped_rows: list[dict[str, Any]] = []
    detailed_runs: list[dict[str, Any]] = []
    for window in windows:
        window_id = window["window_id"]
        sample_ids, variant_ids, positions, _target_alleles, carrier_matrix = load_raw_carriers(
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
                beacon_matrix_full = matrix_for_samples(sample_ids, carrier_matrix, split["beacon"])
                public_ref_matrix_full = matrix_for_samples(sample_ids, carrier_matrix, split["public_ref"])
                template_matrix_full = matrix_for_samples(sample_ids, carrier_matrix, template_ref_ids)
                true_counts_full = beacon_matrix_full.sum(axis=0, dtype=np.int64)

                for profile_index, profile in enumerate(profiles):
                    profile_seed = seeds["group_selection_seed"] + 1009 * profile_index
                    released_indices = _profile_indices(
                        profile,
                        len(variant_ids),
                        beacon_size=beacon_size,
                        seed=profile_seed,
                    )
                    if len(released_indices) == 0:
                        skipped_rows.append(
                            _skipped_row(
                                beacon_size=beacon_size,
                                window_id=window_id,
                                replicate_id=replicate_id,
                                profile=profile,
                                reason=_empty_release_reason(profile),
                                original_variant_count=len(variant_ids),
                            )
                        )
                        utility_rows.extend(
                            _variant_utility_rows(
                                beacon_size=beacon_size,
                                window_id=window_id,
                                replicate_id=replicate_id,
                                profile=profile,
                                variant_ids=variant_ids,
                                true_counts_full=true_counts_full,
                                released_indices=released_indices,
                                released_counts=np.array([], dtype=np.int64),
                            )
                        )
                        continue

                    true_counts = true_counts_full[released_indices]
                    released_counts = _apply_profile(
                        true_counts,
                        profile=profile,
                        beacon_size=beacon_size,
                        seed=profile_seed,
                    )
                    utility_rows.extend(
                        _variant_utility_rows(
                            beacon_size=beacon_size,
                            window_id=window_id,
                            replicate_id=replicate_id,
                            profile=profile,
                            variant_ids=variant_ids,
                            true_counts_full=true_counts_full,
                            released_indices=released_indices,
                            released_counts=released_counts,
                        )
                    )

                    beacon_matrix = beacon_matrix_full[:, released_indices]
                    public_ref_matrix = public_ref_matrix_full[:, released_indices]
                    template_matrix = template_matrix_full[:, released_indices]
                    reference_pairwise = pairwise_moments(public_ref_matrix)
                    split_id = f"{window_id}_replicate_{replicate_id}_N{beacon_size}_{profile}"
                    runtime_config = {
                        "window_id": window_id,
                        "replicate_id": replicate_id,
                        "beacon_size": beacon_size,
                        "release_profile": profile,
                        "config": config,
                    }
                    problem = BeaconAttackProblem(
                        manifest_id=MANIFEST_ID,
                        split_id=split_id,
                        beacon_size=beacon_size,
                        variant_ids=tuple(variant_ids[index] for index in released_indices),
                        positions=tuple(int(positions[index]) for index in released_indices),
                        query_responses=released_counts,
                        target_counts=released_counts,
                        public_ld=reference_pairwise,
                        config=ExperimentConfig(
                            initialization_seed=seeds["initialization_seed"] + 17 * profile_index,
                            optimizer_seed=seeds["optimizer_seed"] + 17 * profile_index,
                            group_selection_seed=seeds["group_selection_seed"],
                            metadata={"window_id": window_id, "experiment": "count_release_perturbation", "profile": profile},
                        ),
                    )
                    truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(split["beacon"]))
                    stats = _release_stats(
                        released_counts,
                        true_counts,
                        released_variant_count=len(released_indices),
                        original_variant_count=len(variant_ids),
                    )

                    method_run_specs: list[dict[str, Any]] = []
                    if not args.saleem_only:
                        method_run_specs.extend(
                            [
                                {
                                    "method": "independent_count",
                                    "reference_ids": [],
                                    "run": lambda problem=problem, reference_pairwise=reference_pairwise, runtime_config=runtime_config, seeds=seeds, profile_index=profile_index: _with_metadata(
                                        independent_count(problem, reference_pairwise=reference_pairwise),
                                        config=runtime_config,
                                        optimizer_seed=seeds["optimizer_seed"] + 17 * profile_index,
                                    ),
                                },
                                {
                                    "method": "pairwise_refinement",
                                    "reference_ids": split["public_ref"],
                                    "run": lambda problem=problem, reference_pairwise=reference_pairwise, runtime_config=runtime_config, seeds=seeds, profile_index=profile_index: _with_metadata(
                                        pairwise_refinement(
                                            problem,
                                            reference_pairwise,
                                            epochs=args.epochs,
                                            proposals_per_epoch=args.proposals_per_epoch,
                                            early_stopping=args.early_stopping,
                                        ),
                                        config=runtime_config,
                                        optimizer_seed=seeds["optimizer_seed"] + 17 * profile_index,
                                    ),
                                },
                            ]
                        )
                    if args.run_saleem:
                        method_run_specs.append(
                            {
                                "method": SALEEM_RECONSTRUCTION_ATTACK,
                                "reference_ids": split["public_ref"],
                                "run": lambda problem=problem, public_ref_matrix=public_ref_matrix, beacon_matrix=beacon_matrix, runtime_config=runtime_config, seeds=seeds, profile_index=profile_index: _with_metadata(
                                    run_saleem_attack(problem, public_ref_matrix, beacon_matrix, args),
                                    config=runtime_config,
                                    optimizer_seed=seeds["optimizer_seed"] + 17 * profile_index,
                                ),
                            }
                        )
                    if not args.saleem_only:
                        method_run_specs.append(
                            {
                                "method": "reference_haplotype_reconstruction",
                                "reference_ids": template_ref_ids,
                                "run": lambda problem=problem, template_matrix=template_matrix, runtime_config=runtime_config, seeds=seeds, profile_index=profile_index: _with_metadata(
                                    reference_haplotype_reconstruction(problem, template_matrix, metadata=runtime_config),
                                    config=runtime_config,
                                    optimizer_seed=seeds["optimizer_seed"] + 17 * profile_index,
                                ),
                            }
                        )

                    for spec in method_run_specs:
                        result, wall_time_sec = run_timed(spec["run"])
                        result = attach_reference_panel_mode(result, shared=args.shared_reference_panel)
                        runtime_fields = runtime_metadata_fields(result, wall_time_sec)
                        evaluation = evaluate_reconstruction(
                            result,
                            truth,
                            target_counts=released_counts,
                            reference_pairwise=reference_pairwise,
                            source_sets={"beacon": split["beacon"], "reference": spec["reference_ids"]},
                            accepted_manifest_ids=None,
                        )
                        row = _record(
                            beacon_size=beacon_size,
                            window_id=window_id,
                            replicate_id=replicate_id,
                            profile=profile,
                            method=spec["method"],
                            result=result,
                            evaluation=evaluation,
                            release_stats=stats,
                            source_overlap_count=len(set(split["beacon"]).intersection(spec["reference_ids"])),
                            runtime_fields=runtime_fields,
                        )
                        rows.append(row)
                        detailed_runs.append(
                            {
                                "summary": row,
                                "metrics": evaluation.metrics.as_dict(),
                                "invalid_reasons": evaluation.invalid_run.reasons,
                                "result_metadata": result.metadata,
                            }
                        )
                    print(f"completed {split_id}")

    rows = merge_saleem_rows(
        new_rows=rows,
        merge_metrics=args.merge_metrics,
        saleem_only=args.saleem_only,
    )
    if args.merge_existing:
        rows = _merge_existing_metrics(
            rows,
            metrics_path=args.out_dir / "count_release_perturbation_metrics.csv",
        )
    aggregate = _aggregate(rows)
    deltas = _profile_deltas(rows)
    _write_csv(args.out_dir / "count_release_perturbation_metrics.csv", rows)
    _write_csv(args.out_dir / "count_release_perturbation_aggregate.csv", aggregate)
    _write_csv(args.out_dir / "count_release_perturbation_deltas.csv", deltas)
    _write_csv(args.out_dir / "count_release_perturbation_variant_utility.csv", utility_rows)
    if skipped_rows:
        _write_csv(args.out_dir / "count_release_perturbation_skipped.csv", skipped_rows)
    _write_json(
        args.out_dir / "count_release_perturbation_manifest.json",
        {"manifest_id": MANIFEST_ID, "config": config, "detailed_runs": detailed_runs},
    )
    _write_report(
        args.out_dir / "count_release_perturbation_report.md",
        aggregate_rows=aggregate,
        delta_rows=deltas,
        skipped_rows=skipped_rows,
        config=config,
    )
    print(f"wrote {args.out_dir / 'count_release_perturbation_report.md'}")


if __name__ == "__main__":
    main()
