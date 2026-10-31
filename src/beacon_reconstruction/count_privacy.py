"""Differentially private release of integer carrier counts.

Each released count is an integer query of sensitivity 1. The mechanism adds
independent discrete Laplace noise and clips the result to the feasible count
range. Clipping is post-processing. Releasing ``m`` counts at per-count budget
``epsilon_per_count`` is ``m * epsilon_per_count``-differentially private by
basic composition.
"""

from __future__ import annotations

from typing import Any

import numpy as np

ACCOUNTING_PER_QUERY = "per_query"
ACCOUNTING_COMPOSED = "composed"
BUDGET_ACCOUNTINGS = (ACCOUNTING_PER_QUERY, ACCOUNTING_COMPOSED)
DEFAULT_NOMINAL_EPSILONS = (0.1, 1.0, 10.0)


def format_epsilon(epsilon: float) -> str:
    """Format a privacy budget without trailing zeros."""

    return f"{float(epsilon):g}"


def calibrate_epsilon_per_count(
    nominal_epsilon: float,
    variant_count: int,
    accounting: str,
) -> float:
    """Return the discrete-Laplace budget used for one carrier count.

    ``per_query`` spends ``nominal_epsilon`` on every count, so the release of
    ``variant_count`` counts has basic-composition budget
    ``nominal_epsilon * variant_count``. ``composed`` splits ``nominal_epsilon``
    evenly across the counts, so the joint release has budget ``nominal_epsilon``.
    """

    if not np.isfinite(nominal_epsilon) or float(nominal_epsilon) <= 0:
        raise ValueError("nominal_epsilon must be positive and finite.")
    if int(variant_count) < 1:
        raise ValueError("variant_count must be positive.")
    if accounting == ACCOUNTING_PER_QUERY:
        return float(nominal_epsilon)
    if accounting == ACCOUNTING_COMPOSED:
        return float(nominal_epsilon) / int(variant_count)
    raise ValueError(f"Unsupported budget accounting: {accounting}.")


def dp_release_profiles(
    nominal_epsilons: tuple[float, ...] = DEFAULT_NOMINAL_EPSILONS,
    accountings: tuple[str, ...] = BUDGET_ACCOUNTINGS,
) -> list[dict[str, Any]]:
    """Return discrete-Laplace profiles for the requested budget accountings."""

    unknown = [accounting for accounting in accountings if accounting not in BUDGET_ACCOUNTINGS]
    if unknown:
        raise ValueError(f"Unsupported budget accounting: {unknown[0]}.")
    profiles: list[dict[str, Any]] = []
    for accounting in accountings:
        for nominal_epsilon in nominal_epsilons:
            profiles.append(
                {
                    "release_profile": f"dp_{accounting}_eps{format_epsilon(nominal_epsilon)}",
                    "budget_accounting": accounting,
                    "nominal_epsilon": float(nominal_epsilon),
                }
            )
    return profiles


def sample_discrete_laplace(
    epsilon: float,
    size: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample the two-sided geometric mechanism with parameter ``epsilon``.

    For sensitivity ``Delta = 1``,
    ``P(Z = k) = ((1 - alpha) / (1 + alpha)) alpha**|k|`` with
    ``alpha = exp(-epsilon)``. A single integer count query protected by this
    noise is ``epsilon``-differentially private.
    """

    count = int(size)
    if count < 0:
        raise ValueError("size must be non-negative.")
    if not np.isfinite(epsilon) or float(epsilon) <= 0:
        raise ValueError("epsilon must be positive and finite.")
    if count == 0:
        return np.zeros(0, dtype=np.int64)

    alpha = float(np.exp(-float(epsilon)))
    if alpha == 0.0:
        return np.zeros(count, dtype=np.int64)
    if alpha >= 1.0:
        raise ValueError("epsilon is too small to sample discrete Laplace noise in float64.")

    zero_probability = (1.0 - alpha) / (1.0 + alpha)
    is_zero = rng.random(count) < zero_probability
    magnitude = rng.geometric(1.0 - alpha, size=count).astype(np.int64, copy=False)
    sign = rng.integers(0, 2, size=count, dtype=np.int64) * 2 - 1
    noise = sign * magnitude
    noise[is_zero] = 0
    return noise


def release_discrete_laplace_counts(
    counts: np.ndarray,
    *,
    beacon_size: int,
    epsilon_per_count: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Add discrete Laplace noise and clip carrier counts to ``[0, beacon_size]``."""

    values = np.asarray(counts)
    if values.ndim != 1:
        raise ValueError("counts must be one-dimensional.")
    if int(beacon_size) < 0:
        raise ValueError("beacon_size must be non-negative.")
    noise = sample_discrete_laplace(epsilon_per_count, values.size, rng)
    released = np.asarray(values, dtype=np.int64) + noise
    return np.clip(released, 0, int(beacon_size)).astype(np.int64, copy=False)
