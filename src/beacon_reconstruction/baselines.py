"""Baseline reconstruction methods for Beacon carrier-count experiments."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .core import BeaconAttackProblem, ReconstructionResult
from .count_feasibility import assert_count_feasible, initialize_independent_count


INDEPENDENT_COUNT = "independent_count"
PAIRWISE_REFINEMENT = "pairwise_refinement"

__all__ = [
    "INDEPENDENT_COUNT",
    "PAIRWISE_REFINEMENT",
    "independent_count",
    "pairwise_loss",
    "pairwise_moments",
    "pairwise_refinement",
]


def _ensure_attack_problem(problem: Any, function_name: str) -> BeaconAttackProblem:
    if not isinstance(problem, BeaconAttackProblem):
        raise TypeError(f"{function_name} requires a BeaconAttackProblem.")
    return problem


def _as_binary_matrix(binary_matrix: ArrayLike, name: str = "binary_matrix") -> NDArray[np.int8]:
    matrix = np.asarray(binary_matrix)
    if matrix.ndim != 2:
        raise ValueError(f"{name} must be a two-dimensional matrix.")
    if not np.all((matrix == 0) | (matrix == 1)):
        raise ValueError(f"{name} must contain only 0/1 values.")
    return matrix.astype(np.int8, copy=True)


def _as_pairwise_reference(reference_pairwise: ArrayLike, variant_count: int) -> NDArray[np.float64]:
    reference = np.asarray(reference_pairwise, dtype=float)
    if reference.shape != (variant_count, variant_count):
        raise ValueError("reference_pairwise must have shape (variant_count, variant_count).")
    if np.any(~np.isfinite(reference)):
        raise ValueError("reference_pairwise must contain only finite values.")
    return reference


def _resolve_seed(explicit_seed: int | None, configured_seed: int | None) -> int:
    if explicit_seed is not None:
        return int(explicit_seed)
    if configured_seed is not None:
        return int(configured_seed)
    return 0


def _max_count_error(matrix: NDArray[np.int8], target_counts: NDArray[np.integer[Any]]) -> int:
    observed = matrix.sum(axis=0, dtype=np.int64)
    return int(np.max(np.abs(observed - np.asarray(target_counts, dtype=np.int64)), initial=0))


def _result_metadata(
    *,
    method_name: str,
    problem: BeaconAttackProblem,
    prior_id: str,
    seeds: dict[str, int],
    max_count_error: int,
    current_pairwise_loss: float | None,
    loss_trace: list[float],
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    metadata: dict[str, Any] = {
        "method_name": method_name,
        "manifest_id": problem.manifest_id,
        "split_id": problem.split_id,
        "prior_id": prior_id,
        "seeds": dict(seeds),
        "max_count_error": max_count_error,
        "pairwise_loss": current_pairwise_loss,
        "loss_trace": list(loss_trace),
    }
    if extra:
        metadata.update(extra)
    return metadata


def pairwise_moments(binary_matrix: ArrayLike) -> NDArray[np.float64]:
    """Return binary co-carrier moments P[j, k] = mean_i B[i, j] B[i, k]."""

    matrix = _as_binary_matrix(binary_matrix)
    n_rows = matrix.shape[0]
    if n_rows == 0:
        raise ValueError("binary_matrix must contain at least one row.")
    return (matrix.T @ matrix).astype(np.float64) / float(n_rows)


def pairwise_loss(binary_matrix: ArrayLike, reference_pairwise: ArrayLike) -> float:
    """Mean squared pairwise-moment error with the diagonal excluded."""

    matrix = _as_binary_matrix(binary_matrix)
    variant_count = matrix.shape[1]
    reference = _as_pairwise_reference(reference_pairwise, variant_count)
    if variant_count < 2:
        return 0.0
    moments = pairwise_moments(matrix)
    upper = np.triu_indices(variant_count, k=1)
    errors = moments[upper] - reference[upper]
    return float(np.mean(errors * errors))


def independent_count(
    problem: BeaconAttackProblem,
    *,
    prior_id: str = "none",
    reference_pairwise: ArrayLike | None = None,
    initialization_seed: int | None = None,
) -> ReconstructionResult:
    """Generate a count-feasible independent carrier assignment."""

    attack_problem = _ensure_attack_problem(problem, "independent_count")
    seed = _resolve_seed(initialization_seed, attack_problem.config.initialization_seed)
    reconstructed = initialize_independent_count(attack_problem, seed=seed)
    assert_count_feasible(attack_problem.target_counts, attack_problem.beacon_size, reconstructed)

    max_count_error = _max_count_error(reconstructed, attack_problem.target_counts)
    current_pairwise_loss = (
        None if reference_pairwise is None else pairwise_loss(reconstructed, reference_pairwise)
    )
    loss_trace = [] if current_pairwise_loss is None else [current_pairwise_loss]
    seeds = {"initialization_seed": seed}
    metadata = _result_metadata(
        method_name=INDEPENDENT_COUNT,
        problem=attack_problem,
        prior_id=prior_id,
        seeds=seeds,
        max_count_error=max_count_error,
        current_pairwise_loss=current_pairwise_loss,
        loss_trace=loss_trace,
    )
    return ReconstructionResult(
        manifest_id=attack_problem.manifest_id,
        split_id=attack_problem.split_id,
        prior_id=prior_id,
        method_name=INDEPENDENT_COUNT,
        seeds=seeds,
        reconstructed_binary=reconstructed,
        metrics={"max_count_error": max_count_error, "pairwise_loss": current_pairwise_loss},
        metadata=metadata,
    )


def _initial_pairwise_state(
    matrix: NDArray[np.int8],
    reference_pairwise: NDArray[np.float64],
) -> tuple[NDArray[np.float64], float, float]:
    variant_count = matrix.shape[1]
    moments = pairwise_moments(matrix)
    if variant_count < 2:
        return moments, 0.0, 0.0
    upper = np.triu_indices(variant_count, k=1)
    loss_sum = float(np.sum((moments[upper] - reference_pairwise[upper]) ** 2))
    return moments, loss_sum, loss_sum / len(upper[0])


def pairwise_refinement(
    problem: BeaconAttackProblem,
    reference_pairwise: ArrayLike,
    *,
    prior_id: str = "pairwise_reference",
    initialization_seed: int | None = None,
    optimizer_seed: int | None = None,
    epochs: int = 20,
    proposals_per_epoch: int = 5000,
    early_stopping: int = 5,
) -> ReconstructionResult:
    """Greedily refine a count-feasible initialization against pairwise moments."""

    attack_problem = _ensure_attack_problem(problem, "pairwise_refinement")
    if epochs < 0:
        raise ValueError("epochs must be non-negative.")
    if proposals_per_epoch < 0:
        raise ValueError("proposals_per_epoch must be non-negative.")
    if early_stopping < 0:
        raise ValueError("early_stopping must be non-negative.")

    init_seed = _resolve_seed(initialization_seed, attack_problem.config.initialization_seed)
    opt_seed = _resolve_seed(optimizer_seed, attack_problem.config.optimizer_seed)
    rng = np.random.default_rng(opt_seed)

    reconstructed = initialize_independent_count(attack_problem, seed=init_seed)
    assert_count_feasible(attack_problem.target_counts, attack_problem.beacon_size, reconstructed)

    n_rows, variant_count = reconstructed.shape
    reference = _as_pairwise_reference(reference_pairwise, variant_count)
    moments, loss_sum, current_loss = _initial_pairwise_state(reconstructed, reference)
    pair_count = variant_count * (variant_count - 1) // 2
    loss_trace = [current_loss]
    accepted_moves = 0
    epochs_without_improvement = 0

    if n_rows > 0 and pair_count > 0:
        target_counts = np.asarray(attack_problem.target_counts, dtype=np.int64)
        eligible_columns = np.flatnonzero((target_counts > 0) & (target_counts < n_rows))
    else:
        eligible_columns = np.array([], dtype=int)

    for _epoch in range(epochs):
        epoch_accepted = 0
        if len(eligible_columns) > 0:
            for _proposal in range(proposals_per_epoch):
                column = int(rng.choice(eligible_columns))
                carrier_rows = np.flatnonzero(reconstructed[:, column] == 1)
                noncarrier_rows = np.flatnonzero(reconstructed[:, column] == 0)
                if len(carrier_rows) == 0 or len(noncarrier_rows) == 0:
                    continue

                row_out = int(rng.choice(carrier_rows))
                row_in = int(rng.choice(noncarrier_rows))
                row_delta = reconstructed[row_in, :].astype(np.float64) - reconstructed[row_out, :].astype(np.float64)
                moment_delta = row_delta / float(n_rows)
                moment_delta[column] = 0.0

                lower_columns = np.arange(column)
                upper_columns = np.arange(column + 1, variant_count)
                lower_old = moments[lower_columns, column] - reference[lower_columns, column]
                lower_new = lower_old + moment_delta[lower_columns]
                upper_old = moments[column, upper_columns] - reference[column, upper_columns]
                upper_new = upper_old + moment_delta[upper_columns]
                delta_loss_sum = float(
                    np.sum(lower_new * lower_new - lower_old * lower_old)
                    + np.sum(upper_new * upper_new - upper_old * upper_old)
                )

                if delta_loss_sum < -1e-15:
                    reconstructed[row_out, column] = 0
                    reconstructed[row_in, column] = 1
                    moments[column, :] += moment_delta
                    moments[:, column] += moment_delta
                    moments[column, column] = float(target_counts[column]) / float(n_rows)
                    loss_sum += delta_loss_sum
                    current_loss = loss_sum / pair_count
                    accepted_moves += 1
                    epoch_accepted += 1

        loss_trace.append(current_loss)
        if epoch_accepted == 0:
            epochs_without_improvement += 1
        else:
            epochs_without_improvement = 0
        if epochs_without_improvement >= early_stopping:
            break

    assert_count_feasible(attack_problem.target_counts, attack_problem.beacon_size, reconstructed)
    max_count_error = _max_count_error(reconstructed, attack_problem.target_counts)
    # Recompute once from the returned matrix to keep reported loss exact.
    final_loss = pairwise_loss(reconstructed, reference)
    loss_trace[-1] = final_loss

    seeds = {"initialization_seed": init_seed, "optimizer_seed": opt_seed}
    metadata = _result_metadata(
        method_name=PAIRWISE_REFINEMENT,
        problem=attack_problem,
        prior_id=prior_id,
        seeds=seeds,
        max_count_error=max_count_error,
        current_pairwise_loss=final_loss,
        loss_trace=loss_trace,
        extra={
            "epochs": int(epochs),
            "proposals_per_epoch": int(proposals_per_epoch),
            "early_stopping": int(early_stopping),
            "accepted_moves": int(accepted_moves),
        },
    )
    return ReconstructionResult(
        manifest_id=attack_problem.manifest_id,
        split_id=attack_problem.split_id,
        prior_id=prior_id,
        method_name=PAIRWISE_REFINEMENT,
        seeds=seeds,
        reconstructed_binary=reconstructed,
        metrics={"max_count_error": max_count_error, "pairwise_loss": final_loss},
        metadata=metadata,
    )
