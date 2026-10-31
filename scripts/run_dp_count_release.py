"""Reconstruction risk under discrete-Laplace carrier-count release.

Each profile adds independent discrete Laplace noise of sensitivity 1 to the
true carrier counts and clips the result to [0, N]. The attack treats those
clipped integers as the released counts. Carrier-level F1 is scored against
the true Beacon matrix.

Two budget accountings are reported for each nominal epsilon:

- per_query spends epsilon on every SNP. Basic composition over m SNPs gives
  a joint budget of m * epsilon.
- composed splits epsilon evenly across the m SNPs, so the joint release has
  budget epsilon and each count uses epsilon / m.
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
from beacon_reconstruction.count_privacy import (
    calibrate_epsilon_per_count,
    dp_release_profiles,
    format_epsilon,
    release_discrete_laplace_counts,
)
from beacon_reconstruction.evaluation import EvaluationResult, evaluate_reconstruction
from beacon_reconstruction.ld_af_baseline import SALEEM_RECONSTRUCTION_ATTACK
from scripts.experiment_runtime_helpers import attach_reference_panel_mode, run_timed, runtime_metadata_fields
from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers, matrix_for_samples as _matrix_for_samples
from scripts.lib.sample_splits import load_sample_manifest, load_seed_manifest, split_samples
from scripts.ld_af_helpers import add_saleem_arguments, run_saleem_attack, saleem_config_dict


MANIFEST_ID = "1000g_eur_hg38_dp_count_release_v1"
DEFAULT_MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
DEFAULT_OUT_DIR = ROOT / "results" / "block_inheritance_v1" / "dp_count_release_shared_panel"
# Fixed indexes preserve the noise seed used by the original full grid.
# New per-query budgets get unused indexes so they do not reuse an old draw.
NOISE_PROFILE_INDEX = {
    "exact_count": 0,
    "dp_per_query_eps0.1": 1,
    "dp_per_query_eps1": 2,
    "dp_per_query_eps10": 3,
    "dp_composed_eps0.1": 4,
    "dp_composed_eps1": 5,
    "dp_composed_eps10": 6,
    "dp_per_query_eps0.01": 7,
}
DEFAULT_METHODS = (
    "independent_count",
    "pairwise_refinement",
    "reference_haplotype_reconstruction",
    SALEEM_RECONSTRUCTION_ATTACK,
)
_FLOAT_FIELDS = (
    "precision",
    "recall",
    "f1",
    "alignment_total_hamming_cost",
    "max_count_error_to_released_counts",
    "count_mae_to_true",
    "changed_count_fraction",
    "boundary_count_fraction",
    "wall_time_sec",
    "source_overlap_count",
    "nominal_epsilon",
    "epsilon_per_count",
    "composed_epsilon",
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


def _release_stats(
    released_counts: np.ndarray,
    true_counts: np.ndarray,
    *,
    beacon_size: int,
) -> dict[str, float | int]:
    released = np.asarray(released_counts, dtype=np.int64)
    delta = released - np.asarray(true_counts, dtype=np.int64)
    on_boundary = (released == 0) | (released == int(beacon_size))
    return {
        "count_mae_to_true": float(np.mean(np.abs(delta))) if delta.size else float("nan"),
        "count_max_abs_error": int(np.max(np.abs(delta), initial=0)),
        "changed_count_fraction": float(np.mean(delta != 0)) if delta.size else float("nan"),
        "boundary_count_fraction": float(np.mean(on_boundary)) if on_boundary.size else float("nan"),
    }


def _mean(values: list[float]) -> float:
    return float(np.mean(values)) if values else float("nan")


def _std(values: list[float]) -> float:
    return float(np.std(values, ddof=1)) if len(values) > 1 else float("nan")


def _aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = sorted(
        {
            (
                int(row["beacon_size"]),
                row["release_profile"],
                row["method"],
            )
            for row in rows
        },
        key=lambda item: (item[0], _profile_order(item[1]), item[2]),
    )
    out: list[dict[str, Any]] = []
    for beacon_size, profile, method in keys:
        selected = [
            row
            for row in rows
            if int(row["beacon_size"]) == beacon_size and row["release_profile"] == profile and row["method"] == method
        ]
        f1_values = [float(row["f1"]) for row in selected]
        out.append(
            {
                "beacon_size": beacon_size,
                "release_profile": profile,
                "budget_accounting": selected[0]["budget_accounting"],
                "nominal_epsilon": selected[0]["nominal_epsilon"],
                "epsilon_per_count": selected[0]["epsilon_per_count"],
                "composed_epsilon": selected[0]["composed_epsilon"],
                "method": method,
                "n_runs": len(selected),
                "mean_f1": _mean(f1_values),
                "sd_f1": _std(f1_values),
                "mean_count_mae_to_true": _mean([float(row["count_mae_to_true"]) for row in selected]),
                "mean_changed_count_fraction": _mean([float(row["changed_count_fraction"]) for row in selected]),
                "mean_boundary_count_fraction": _mean([float(row["boundary_count_fraction"]) for row in selected]),
                "mean_max_count_error_to_released_counts": _mean(
                    [float(row["max_count_error_to_released_counts"]) for row in selected]
                ),
                "invalid_run_count": int(sum(bool(row["invalid_run"]) for row in selected)),
                "mean_wall_time_sec": _mean([float(row["wall_time_sec"]) for row in selected]),
            }
        )
    return out


def _profile_order(profile: str) -> tuple[int, float]:
    if profile == "exact_count":
        return (0, 0.0)
    if profile.startswith("dp_per_query_eps"):
        return (1, float(profile.removeprefix("dp_per_query_eps")))
    if profile.startswith("dp_composed_eps"):
        return (2, float(profile.removeprefix("dp_composed_eps")))
    return (3, 0.0)


def _profile_deltas(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {
        (int(row["beacon_size"]), row["window_id"], row["replicate_id"], row["release_profile"], row["method"]): row
        for row in rows
    }
    units = sorted({(int(row["beacon_size"]), row["window_id"], row["replicate_id"]) for row in rows})
    methods = sorted({row["method"] for row in rows})
    profiles = sorted({row["release_profile"] for row in rows} - {"exact_count"}, key=_profile_order)
    out: list[dict[str, Any]] = []
    for beacon_size, window_id, replicate_id in units:
        for profile in profiles:
            for method in methods:
                exact_key = (beacon_size, window_id, replicate_id, "exact_count", method)
                profile_key = (beacon_size, window_id, replicate_id, profile, method)
                if exact_key not in by_key or profile_key not in by_key:
                    continue
                exact = by_key[exact_key]
                private = by_key[profile_key]
                out.append(
                    {
                        "beacon_size": beacon_size,
                        "window_id": window_id,
                        "replicate_id": replicate_id,
                        "release_profile": profile,
                        "budget_accounting": private["budget_accounting"],
                        "nominal_epsilon": private["nominal_epsilon"],
                        "method": method,
                        "delta_f1_vs_exact_count": float(private["f1"]) - float(exact["f1"]),
                    }
                )
    return out


def _aggregate_deltas(delta_rows: list[dict[str, Any]]) -> dict[tuple[int, str, str], float]:
    grouped: dict[tuple[int, str, str], list[float]] = {}
    for row in delta_rows:
        key = (int(row["beacon_size"]), row["release_profile"], row["method"])
        grouped.setdefault(key, []).append(float(row["delta_f1_vs_exact_count"]))
    return {key: _mean(values) for key, values in grouped.items()}


def _write_report(
    path: Path,
    *,
    aggregate_rows: list[dict[str, Any]],
    delta_means: dict[tuple[int, str, str], float],
    config: dict[str, Any],
) -> None:
    lines = [
        "# Discrete-Laplace Count Release",
        "",
        "Carrier counts are released with the discrete Laplace mechanism at sensitivity 1.",
        "Noise is independent across SNPs and the released integer is clipped to [0, N].",
        "Clipping is post-processing. The joint privacy claim uses basic composition,",
        "not advanced composition and not a correlation-aware accountant.",
        "",
        "Per-query accounting spends the nominal epsilon on every SNP, so m SNPs compose to m * epsilon.",
        "Composed accounting splits the nominal epsilon evenly, so each SNP uses epsilon / m and the joint budget is epsilon.",
        "Reconstruction is scored against the true carrier matrix. The attack only observes the clipped noisy counts.",
        "",
        "## Configuration",
        "",
        f"- windows: {config['max_windows']}",
        f"- replicates: {config['max_replicates']}",
        f"- beacon sizes: {', '.join(str(size) for size in config['beacon_sizes'])}",
        f"- nominal epsilons: {', '.join(format_epsilon(value) for value in config['nominal_epsilons'])}",
        f"- methods: {', '.join(config['methods'])}",
        f"- reference panel: {config['reference_panel_mode']}",
        f"- template selection: {config['template_selection_rule']}",
        f"- pairwise epochs: {config['pairwise_epochs']}",
        f"- pairwise proposals per epoch: {config['pairwise_proposals_per_epoch']}",
        f"- Saleem correlation/frequency epochs: {config['saleem_corr_epoch']}/{config['saleem_freq_epoch']}",
        f"- Saleem max outer iterations: {config['saleem_max_iterations']}",
        "",
        "## Aggregate",
        "",
        "| N | profile | accounting | nominal ε | ε per count | composed ε | method | runs | F1 | SD | ΔF1 | count MAE | changed |",
        "|---:|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        delta = delta_means.get((int(row["beacon_size"]), row["release_profile"], row["method"]))
        delta_text = "" if delta is None or math.isnan(delta) else f"{delta:.4f}"
        nominal = row["nominal_epsilon"]
        per_count = row["epsilon_per_count"]
        composed = row["composed_epsilon"]
        lines.append(
            f"| {row['beacon_size']} | {row['release_profile']} | {row['budget_accounting']} | "
            f"{'' if nominal == '' else format_epsilon(float(nominal))} | "
            f"{'' if per_count == '' else f'{float(per_count):.6g}'} | "
            f"{'' if composed == '' else format_epsilon(float(composed))} | "
            f"{row['method']} | {row['n_runs']} | {row['mean_f1']:.4f} | {row['sd_f1']:.4f} | "
            f"{delta_text} | {row['mean_count_mae_to_true']:.4f} | {row['mean_changed_count_fraction']:.4f} |"
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Per-query epsilon is the privacy loss of one count. It is not the privacy loss of the 500-SNP release.",
            "Composed epsilon is the basic-composition budget of the whole count vector.",
            "Rows that match the noisy counts can still differ from the true genotypes; F1 measures that difference.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def _row_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        int(row["beacon_size"]),
        str(row["window_id"]),
        str(row["replicate_id"]),
        str(row["release_profile"]),
        str(row["method"]),
    )


def _coerce_row(row: dict[str, Any]) -> dict[str, Any]:
    out = dict(row)
    out["beacon_size"] = int(row["beacon_size"])
    for key in _FLOAT_FIELDS:
        value = row.get(key, "")
        if value == "" or value is None:
            continue
        out[key] = float(value)
    invalid = row.get("invalid_run", False)
    if isinstance(invalid, str):
        out["invalid_run"] = invalid.strip().lower() in {"1", "true", "yes"}
    return out


def _reference_ids(method: str, public_ref: list[str], template_ref: list[str]) -> list[str]:
    if method == "independent_count":
        return []
    if method == "reference_haplotype_reconstruction":
        return template_ref
    return public_ref


def _run_method(
    method: str,
    *,
    problem: BeaconAttackProblem,
    reference_pairwise: np.ndarray,
    template_matrix: np.ndarray,
    public_ref_matrix: np.ndarray,
    beacon_matrix: np.ndarray,
    runtime_config: dict[str, Any],
    optimizer_seed: int,
    pairwise_epochs: int,
    pairwise_proposals_per_epoch: int,
    pairwise_early_stopping: int,
    saleem_args: argparse.Namespace,
) -> ReconstructionResult:
    if method == "independent_count":
        result = independent_count(problem, reference_pairwise=reference_pairwise)
    elif method == "pairwise_refinement":
        result = pairwise_refinement(
            problem,
            reference_pairwise,
            epochs=pairwise_epochs,
            proposals_per_epoch=pairwise_proposals_per_epoch,
            early_stopping=pairwise_early_stopping,
        )
    elif method == "reference_haplotype_reconstruction":
        result = reference_haplotype_reconstruction(problem, template_matrix, metadata=runtime_config)
    elif method == SALEEM_RECONSTRUCTION_ATTACK:
        result = run_saleem_attack(problem, public_ref_matrix, beacon_matrix, saleem_args)
    else:
        raise ValueError(f"Unsupported method: {method}.")
    return _with_metadata(result, config=runtime_config, optimizer_seed=optimizer_seed)


def _publish(out_dir: Path, rows: list[dict[str, Any]], config: dict[str, Any]) -> None:
    aggregate = _aggregate(rows)
    deltas = _profile_deltas(rows)
    delta_means = _aggregate_deltas(deltas)
    _write_csv(out_dir / "dp_count_release_metrics.csv", rows)
    _write_csv(out_dir / "dp_count_release_aggregate.csv", aggregate)
    if deltas:
        _write_csv(out_dir / "dp_count_release_deltas.csv", deltas)
    _write_json(out_dir / "dp_count_release_manifest.json", {"manifest_id": MANIFEST_ID, "config": config})
    _write_report(
        out_dir / "dp_count_release_report.md",
        aggregate_rows=aggregate,
        delta_means=delta_means,
        config=config,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--beacon-sizes", default="25,50,100,250")
    parser.add_argument("--epsilons", default="0.01,0.1,1")
    parser.add_argument("--accountings", default="per_query", help="Comma-separated: per_query, composed.")
    parser.add_argument("--skip-exact", action="store_true", help="Do not run the unmodified count release.")
    parser.add_argument("--methods", default=",".join(DEFAULT_METHODS))
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--proposals-per-epoch", type=int, default=5000)
    parser.add_argument("--early-stopping", type=int, default=5)
    parser.add_argument("--multiwindow-dir", type=Path, default=DEFAULT_MULTIWINDOW_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Keep completed rows in the metrics CSV and run only missing method profiles.",
    )
    add_saleem_arguments(parser)
    args = parser.parse_args()

    nominal_epsilons = tuple(float(value) for value in args.epsilons.split(",") if value.strip())
    if any(value <= 0 for value in nominal_epsilons):
        parser.error("epsilons must be positive")
    methods = [value.strip() for value in args.methods.split(",") if value.strip()]
    unknown = sorted(set(methods) - set(DEFAULT_METHODS))
    if unknown:
        parser.error(f"unsupported methods: {', '.join(unknown)}")
    beacon_sizes = sorted({int(value) for value in args.beacon_sizes.split(",") if value.strip()})
    accountings = tuple(value.strip() for value in args.accountings.split(",") if value.strip())
    try:
        profiles = dp_release_profiles(nominal_epsilons, accountings=accountings)
    except ValueError as exc:
        parser.error(str(exc))
    release_plan: list[dict[str, Any]] = []
    if not args.skip_exact:
        release_plan.append({"release_profile": "exact_count", "budget_accounting": "", "nominal_epsilon": ""})
    release_plan.extend(profiles)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    windows = _read_csv(args.multiwindow_dir / "multiwindow_manifest.csv")[: max(0, args.max_windows)]
    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]
    config = {
        "max_windows": args.max_windows,
        "max_replicates": args.max_replicates,
        "beacon_sizes": beacon_sizes,
        "nominal_epsilons": list(nominal_epsilons),
        "release_profiles": [item["release_profile"] for item in release_plan],
        "methods": methods,
        "multiwindow_dir": str(args.multiwindow_dir),
        "out_dir": str(args.out_dir),
        "reference_panel_mode": "shared",
        "template_selection_rule": TEMPLATE_SELECTION_RULE,
        "mechanism": "discrete_laplace",
        "sensitivity": 1,
        "composition": "basic",
        "pairwise_epochs": args.epochs,
        "pairwise_proposals_per_epoch": args.proposals_per_epoch,
        "pairwise_early_stopping": args.early_stopping,
        **saleem_config_dict(args),
    }

    metrics_path = args.out_dir / "dp_count_release_metrics.csv"
    rows = [_coerce_row(row) for row in _read_csv(metrics_path)] if args.resume and metrics_path.exists() else []
    completed = {_row_key(row) for row in rows}
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
                    shared_reference_panel=True,
                )
                template_ref_ids = split["public_ref"]
                beacon_matrix = _matrix_for_samples(sample_ids, carrier_matrix, split["beacon"])
                public_ref_matrix = _matrix_for_samples(sample_ids, carrier_matrix, split["public_ref"])
                template_matrix = _matrix_for_samples(sample_ids, carrier_matrix, template_ref_ids)
                true_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                reference_pairwise = pairwise_moments(public_ref_matrix)
                variant_count = int(len(variant_ids))

                for profile in release_plan:
                    profile_name = str(profile["release_profile"])
                    accounting = str(profile["budget_accounting"])
                    nominal_epsilon = profile["nominal_epsilon"]
                    if profile_name == "exact_count":
                        released_counts = true_counts.copy()
                        epsilon_per_count = ""
                        composed_epsilon = ""
                    else:
                        epsilon_value = calibrate_epsilon_per_count(
                            float(nominal_epsilon),
                            variant_count,
                            accounting,
                        )
                        if profile_name not in NOISE_PROFILE_INDEX:
                            raise ValueError(f"No stable noise index for profile {profile_name}.")
                        noise_seed = int(seeds["group_selection_seed"]) + 1009 * NOISE_PROFILE_INDEX[profile_name]
                        released_counts = release_discrete_laplace_counts(
                            true_counts,
                            beacon_size=beacon_size,
                            epsilon_per_count=epsilon_value,
                            rng=np.random.default_rng(noise_seed),
                        )
                        epsilon_per_count = epsilon_value
                        composed_epsilon = epsilon_value * variant_count
                    split_id = f"{window_id}_replicate_{replicate_id}_N{beacon_size}_{profile_name}"
                    pending_methods = [
                        method
                        for method in methods
                        if (beacon_size, window_id, str(replicate_id), profile_name, method) not in completed
                    ]
                    if not pending_methods:
                        print(f"skipped {split_id}", flush=True)
                        continue
                    stats = _release_stats(released_counts, true_counts, beacon_size=beacon_size)
                    optimizer_seed = int(seeds["optimizer_seed"])
                    runtime_config = {
                        "window_id": window_id,
                        "replicate_id": replicate_id,
                        "beacon_size": beacon_size,
                        "release_profile": profile_name,
                        "budget_accounting": accounting,
                        "nominal_epsilon": nominal_epsilon,
                        "epsilon_per_count": epsilon_per_count,
                        "composed_epsilon": composed_epsilon,
                        "config": config,
                    }
                    problem = BeaconAttackProblem(
                        manifest_id=MANIFEST_ID,
                        split_id=split_id,
                        beacon_size=beacon_size,
                        variant_ids=tuple(variant_ids),
                        positions=tuple(int(position) for position in positions),
                        query_responses=released_counts,
                        target_counts=released_counts,
                        public_ld=reference_pairwise,
                        config=ExperimentConfig(
                            initialization_seed=int(seeds["initialization_seed"]),
                            optimizer_seed=optimizer_seed,
                            group_selection_seed=int(seeds["group_selection_seed"]),
                            metadata={
                                "window_id": window_id,
                                "experiment": "dp_count_release",
                                "profile": profile_name,
                            },
                        ),
                    )
                    truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(split["beacon"]))
                    for method in pending_methods:
                        result, wall_time_sec = run_timed(
                            lambda method=method, problem=problem, reference_pairwise=reference_pairwise, template_matrix=template_matrix, public_ref_matrix=public_ref_matrix, beacon_matrix=beacon_matrix, runtime_config=runtime_config, optimizer_seed=optimizer_seed: _run_method(
                                method,
                                problem=problem,
                                reference_pairwise=reference_pairwise,
                                template_matrix=template_matrix,
                                public_ref_matrix=public_ref_matrix,
                                beacon_matrix=beacon_matrix,
                                runtime_config=runtime_config,
                                optimizer_seed=optimizer_seed,
                                pairwise_epochs=args.epochs,
                                pairwise_proposals_per_epoch=args.proposals_per_epoch,
                                pairwise_early_stopping=args.early_stopping,
                                saleem_args=args,
                            )
                        )
                        result = attach_reference_panel_mode(result, shared=True)
                        reference_ids = _reference_ids(method, split["public_ref"], template_ref_ids)
                        evaluation: EvaluationResult = evaluate_reconstruction(
                            result,
                            truth,
                            target_counts=released_counts,
                            reference_pairwise=reference_pairwise,
                            source_sets={"beacon": split["beacon"], "reference": reference_ids},
                            accepted_manifest_ids=None,
                        )
                        metrics = evaluation.metrics
                        rows.append(
                            {
                                "beacon_size": beacon_size,
                                "window_id": window_id,
                                "replicate_id": replicate_id,
                                "release_profile": profile_name,
                                "budget_accounting": accounting,
                                "nominal_epsilon": nominal_epsilon,
                                "epsilon_per_count": epsilon_per_count,
                                "composed_epsilon": composed_epsilon,
                                "variant_count": variant_count,
                                "method": method,
                                "method_name": result.method_name,
                                "precision": metrics.precision,
                                "recall": metrics.recall,
                                "f1": metrics.f1,
                                "alignment_total_hamming_cost": metrics.alignment_total_hamming_cost,
                                "max_count_error_to_released_counts": metrics.max_count_error,
                                "source_overlap_count": len(set(split["beacon"]).intersection(reference_ids)),
                                **stats,
                                "invalid_run": evaluation.invalid_run.is_invalid,
                                "invalid_reasons": ";".join(evaluation.invalid_run.reasons),
                                **runtime_metadata_fields(result, wall_time_sec),
                            }
                        )
                        completed.add(_row_key(rows[-1]))
                        _publish(args.out_dir, rows, config)
                        print(
                            f"completed {split_id} {method} f1={metrics.f1:.4f} wall={wall_time_sec:.1f}s",
                            flush=True,
                        )
                    print(f"completed {split_id}", flush=True)

    _publish(args.out_dir, rows, config)
    print(f"wrote {args.out_dir / 'dp_count_release_report.md'}")


if __name__ == "__main__":
    main()
