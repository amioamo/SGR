from __future__ import annotations

from itertools import combinations_with_replacement

import numpy as np

from beacon_reconstruction.reconstruction import _select_reference_templates, minimum_template_count_l1


def _brute_force_l1(reference: np.ndarray, target: np.ndarray, beacon_size: int) -> int:
    best = None
    for combo in combinations_with_replacement(range(reference.shape[0]), beacon_size):
        column_sums = reference[list(combo)].sum(axis=0, dtype=np.int64)
        value = int(np.abs(column_sums - target).sum())
        best = value if best is None else min(best, value)
    if best is None:
        raise AssertionError("empty enumeration")
    return best


def test_exact_minimum_matches_independent_enumeration() -> None:
    reference = np.array(
        [
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
        ],
        dtype=np.int8,
    )
    target = np.array([2, 2, 2], dtype=np.int64)

    assert minimum_template_count_l1(reference, target, beacon_size=3) == _brute_force_l1(reference, target, 3)
    assert minimum_template_count_l1(reference, target, beacon_size=3) == 0


def test_greedy_gap_is_nonnegative_and_can_miss_the_optimum() -> None:
    reference = np.array(
        [
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
        ],
        dtype=np.int8,
    )
    target = np.array([2, 2, 2], dtype=np.int64)
    greedy = _select_reference_templates(reference, target, beacon_size=3, rng=np.random.default_rng(0))
    greedy_l1 = int(np.abs(greedy.sum(axis=0, dtype=np.int64) - target).sum())
    optimal_l1 = minimum_template_count_l1(reference, target, beacon_size=3)

    assert optimal_l1 == 0
    assert greedy_l1 >= optimal_l1
    assert greedy_l1 - optimal_l1 == 1
    assert (greedy_l1 - optimal_l1) / max(1, optimal_l1) == 1


def test_greedy_hits_the_optimum_when_every_reference_row_is_identical() -> None:
    reference = np.ones((4, 3), dtype=np.int8)
    target = np.array([2, 2, 2], dtype=np.int64)
    greedy = _select_reference_templates(reference, target, beacon_size=2, rng=np.random.default_rng(1))
    greedy_l1 = int(np.abs(greedy.sum(axis=0, dtype=np.int64) - target).sum())
    optimal_l1 = minimum_template_count_l1(reference, target, beacon_size=2)

    assert greedy_l1 == optimal_l1 == 0
