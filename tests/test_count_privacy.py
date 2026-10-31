from __future__ import annotations

import numpy as np
import pytest

from beacon_reconstruction.count_privacy import (
    ACCOUNTING_COMPOSED,
    ACCOUNTING_PER_QUERY,
    calibrate_epsilon_per_count,
    dp_release_profiles,
    release_discrete_laplace_counts,
    sample_discrete_laplace,
)


def test_budget_calibration_splits_only_under_composition() -> None:
    assert calibrate_epsilon_per_count(1.0, 500, ACCOUNTING_PER_QUERY) == 1.0
    assert calibrate_epsilon_per_count(1.0, 500, ACCOUNTING_COMPOSED) == pytest.approx(0.002)
    assert calibrate_epsilon_per_count(0.1, 500, ACCOUNTING_COMPOSED) == pytest.approx(0.0002)


def test_budget_calibration_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError):
        calibrate_epsilon_per_count(0.0, 500, ACCOUNTING_PER_QUERY)
    with pytest.raises(ValueError):
        calibrate_epsilon_per_count(1.0, 0, ACCOUNTING_PER_QUERY)
    with pytest.raises(ValueError):
        calibrate_epsilon_per_count(1.0, 500, "advanced")


def test_per_query_profiles_can_use_a_finer_epsilon_grid() -> None:
    profiles = dp_release_profiles((0.01, 0.1, 1.0), accountings=("per_query",))
    assert [profile["release_profile"] for profile in profiles] == [
        "dp_per_query_eps0.01",
        "dp_per_query_eps0.1",
        "dp_per_query_eps1",
    ]


def test_profile_grid_covers_both_accountings() -> None:
    profiles = dp_release_profiles()
    assert [profile["release_profile"] for profile in profiles] == [
        "dp_per_query_eps0.1",
        "dp_per_query_eps1",
        "dp_per_query_eps10",
        "dp_composed_eps0.1",
        "dp_composed_eps1",
        "dp_composed_eps10",
    ]


def test_discrete_laplace_matches_two_sided_geometric_probabilities() -> None:
    epsilon = float(np.log(2.0))
    alpha = 0.5
    scale = (1.0 - alpha) / (1.0 + alpha)
    samples = sample_discrete_laplace(epsilon, 200_000, np.random.default_rng(7))
    for value in range(-3, 4):
        expected = scale * alpha ** abs(value)
        observed = float(np.mean(samples == value))
        assert observed == pytest.approx(expected, abs=0.004)


def test_discrete_laplace_is_deterministic_given_the_seed() -> None:
    first = sample_discrete_laplace(1.0, 30, np.random.default_rng(11))
    second = sample_discrete_laplace(1.0, 30, np.random.default_rng(11))
    assert np.array_equal(first, second)


def test_released_counts_stay_inside_the_cohort_bounds() -> None:
    rng = np.random.default_rng(3)
    counts = np.array([0, 0, 4, 10, 10], dtype=np.int64)
    released = release_discrete_laplace_counts(
        np.tile(counts, 200),
        beacon_size=10,
        epsilon_per_count=0.05,
        rng=rng,
    )
    assert released.dtype == np.int64
    assert int(released.min()) >= 0
    assert int(released.max()) <= 10


def test_large_epsilon_adds_almost_no_noise() -> None:
    counts = np.arange(20, dtype=np.int64)
    released = release_discrete_laplace_counts(
        counts,
        beacon_size=50,
        epsilon_per_count=20.0,
        rng=np.random.default_rng(5),
    )
    assert np.array_equal(released, counts)
