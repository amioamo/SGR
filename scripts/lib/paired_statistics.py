"""Paired delta summaries used by experiment scripts."""

from __future__ import annotations

import math
from typing import Any

import numpy as np


def _two_sided_sign_test_p(wins: int, losses: int) -> float:
    n = wins + losses
    if n == 0:
        return float("nan")
    smaller = min(wins, losses)
    tail = sum(math.comb(n, k) for k in range(smaller + 1)) / (2**n)
    return min(1.0, 2.0 * tail)


def _wilcoxon_p(deltas: np.ndarray) -> float:
    nonzero = deltas[np.abs(deltas) > 1e-12]
    if nonzero.size == 0:
        return float("nan")
    try:
        from scipy.stats import wilcoxon
    except ImportError:
        return float("nan")
    try:
        return float(wilcoxon(nonzero, alternative="two-sided").pvalue)
    except ValueError:
        return float("nan")


def _bootstrap_mean_ci(deltas: np.ndarray, *, seed: int, iterations: int) -> tuple[float, float]:
    finite = deltas[np.isfinite(deltas)]
    if finite.size == 0:
        return float("nan"), float("nan")
    if finite.size == 1 or iterations <= 0:
        value = float(np.mean(finite))
        return value, value
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, finite.size, size=(iterations, finite.size))
    means = finite[indices].mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def _summarize_deltas(
    *,
    comparison: str,
    metric: str,
    deltas: list[float],
    baseline_values: list[float],
    method_values: list[float],
    higher_is_better: bool,
    bootstrap_seed: int,
    bootstrap_iterations: int,
) -> dict[str, Any]:
    delta_array = np.asarray(deltas, dtype=float)
    baseline_array = np.asarray(baseline_values, dtype=float)
    method_array = np.asarray(method_values, dtype=float)
    finite_delta = delta_array[np.isfinite(delta_array)]
    if higher_is_better:
        wins = int(np.sum(finite_delta > 1e-12))
        ties = int(np.sum(np.abs(finite_delta) <= 1e-12))
        losses = int(np.sum(finite_delta < -1e-12))
    else:
        wins = int(np.sum(finite_delta < -1e-12))
        ties = int(np.sum(np.abs(finite_delta) <= 1e-12))
        losses = int(np.sum(finite_delta > 1e-12))
    ci_low, ci_high = _bootstrap_mean_ci(
        finite_delta,
        seed=bootstrap_seed,
        iterations=bootstrap_iterations,
    )
    return {
        "comparison": comparison,
        "metric": metric,
        "n_pairs": int(finite_delta.size),
        "mean_baseline": float(np.nanmean(baseline_array)) if baseline_array.size else float("nan"),
        "mean_method": float(np.nanmean(method_array)) if method_array.size else float("nan"),
        "mean_delta": float(np.nanmean(delta_array)) if delta_array.size else float("nan"),
        "median_delta": float(np.nanmedian(delta_array)) if delta_array.size else float("nan"),
        "ci95_low": ci_low,
        "ci95_high": ci_high,
        "wins": wins,
        "ties": ties,
        "losses": losses,
        "sign_test_p": _two_sided_sign_test_p(wins, losses),
        "wilcoxon_p": _wilcoxon_p(finite_delta),
        "higher_is_better": higher_is_better,
    }

