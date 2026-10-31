"""Carrier-count feasibility utilities for reconstruction methods."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .core import BeaconAttackProblem


def _ensure_attack_problem(problem: Any, function_name: str) -> BeaconAttackProblem:
    if not isinstance(problem, BeaconAttackProblem):
        raise TypeError(f"{function_name} requires a BeaconAttackProblem.")
    return problem


def _as_binary_matrix(matrix: ArrayLike, name: str = "carrier_matrix") -> NDArray[np.int8]:
    array = np.asarray(matrix)
    if array.ndim != 2:
        raise ValueError(f"{name} must be a two-dimensional matrix.")
    if not np.all((array == 0) | (array == 1)):
        raise ValueError(f"{name} must contain only 0/1 values.")
    return array.astype(np.int8, copy=True)


def assert_count_feasible(
    target_counts: ArrayLike,
    beacon_size: int,
    carrier_matrix: ArrayLike | None = None,
) -> NDArray[np.int64]:
    """Validate hard carrier counts and, optionally, a matrix against those counts."""

    if beacon_size < 0:
        raise ValueError("beacon_size must be non-negative.")
    counts = np.asarray(target_counts)
    if counts.ndim != 1:
        raise ValueError("target_counts must be one-dimensional.")
    if not np.issubdtype(counts.dtype, np.integer):
        if np.any(~np.isfinite(counts.astype(float))):
            raise ValueError("target_counts must contain finite integer values.")
        if not np.all(np.equal(counts, np.floor(counts))):
            raise ValueError("target_counts must contain integer carrier counts.")
        counts = counts.astype(np.int64)
    else:
        counts = counts.astype(np.int64, copy=True)
    if np.any(counts < 0) or np.any(counts > beacon_size):
        raise ValueError("target_counts must be between 0 and beacon_size.")

    if carrier_matrix is not None:
        matrix = _as_binary_matrix(carrier_matrix)
        expected_shape = (beacon_size, len(counts))
        if matrix.shape != expected_shape:
            raise ValueError(f"carrier_matrix must have shape {expected_shape}.")
        observed = matrix.sum(axis=0, dtype=np.int64)
        if not np.array_equal(observed, counts):
            raise ValueError("carrier_matrix column sums do not match target_counts.")

    return counts


def initialize_independent_count(
    problem: BeaconAttackProblem,
    seed: int | None = None,
) -> NDArray[np.int8]:
    """Sample an independent binary matrix that exactly matches target carrier counts."""

    attack_problem = _ensure_attack_problem(problem, "initialize_independent_count")
    if attack_problem.response_channel != "carrier_count":
        raise ValueError("Only carrier_count responses define hard target counts.")
    target_counts = np.asarray(attack_problem.target_counts, dtype=np.int64)
    assert_count_feasible(target_counts, attack_problem.beacon_size)
    initialization_seed = attack_problem.config.initialization_seed if seed is None else seed
    rng = np.random.default_rng(initialization_seed)

    matrix = np.zeros((attack_problem.beacon_size, len(target_counts)), dtype=np.int8)
    for column_index, carrier_count in enumerate(target_counts):
        if carrier_count == 0:
            continue
        carrier_rows = rng.choice(attack_problem.beacon_size, size=int(carrier_count), replace=False)
        matrix[carrier_rows, column_index] = 1

    assert_count_feasible(target_counts, attack_problem.beacon_size, matrix)
    return matrix
