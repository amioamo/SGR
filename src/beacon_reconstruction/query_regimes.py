"""SNP-panel selectors for LD strength, SNP density, and genomic span.

Every regime returns a fixed number of variants. Density and span are then
properties of how those variants are placed, and LD pruning changes which
variants survive inside a longer candidate pool.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RegimeSpec:
    regime: str
    description: str
    r2_threshold: float | None = None
    target_gap_bp: int | None = None
    target_span_bp: int | None = None


REGIME_SPECS: tuple[RegimeSpec, ...] = (
    RegimeSpec(
        "dense",
        "consecutive candidate SNPs; highest density and shortest span",
    ),
    RegimeSpec(
        "ld_r2_080",
        "genomic-order greedy prune with pairwise carrier r^2 <= 0.80",
        r2_threshold=0.80,
    ),
    RegimeSpec(
        "ld_r2_050",
        "genomic-order greedy prune with pairwise carrier r^2 <= 0.50",
        r2_threshold=0.50,
    ),
    RegimeSpec(
        "ld_r2_010",
        "genomic-order greedy prune with pairwise carrier r^2 <= 0.10",
        r2_threshold=0.10,
    ),
    RegimeSpec(
        "sparse",
        "low SNP density with a minimum inter-SNP gap of 50 kb",
        target_gap_bp=50_000,
    ),
    RegimeSpec(
        "long_range",
        "fixed SNP count evenly spaced across a 10 Mb span",
        target_span_bp=10_000_000,
    ),
)


def regime_by_name(name: str) -> RegimeSpec:
    for spec in REGIME_SPECS:
        if spec.regime == name:
            return spec
    known = ", ".join(spec.regime for spec in REGIME_SPECS)
    raise KeyError(f"Unknown regime {name!r}. Known regimes: {known}.")


def select_consecutive(start_index: int, count: int, n_variants: int) -> np.ndarray:
    if count <= 0:
        raise ValueError("count must be positive.")
    end = start_index + count
    if start_index < 0 or end > n_variants:
        raise ValueError(f"Consecutive window [{start_index}, {end}) does not fit in {n_variants} variants.")
    return np.arange(start_index, end, dtype=int)


def select_target_gap(
    positions: np.ndarray,
    eligible_indices: np.ndarray,
    *,
    count: int,
    target_gap_bp: int,
) -> np.ndarray:
    """Take `count` eligible variants whose successive positions are at least `target_gap_bp` apart."""
    if count <= 0:
        raise ValueError("count must be positive.")
    if target_gap_bp <= 0:
        raise ValueError("target_gap_bp must be positive.")
    if len(eligible_indices) == 0:
        raise ValueError("No eligible variants.")
    selected = [int(eligible_indices[0])]
    last_pos = int(positions[selected[0]])
    for index in eligible_indices[1:]:
        pos = int(positions[int(index)])
        if pos - last_pos >= target_gap_bp:
            selected.append(int(index))
            last_pos = pos
            if len(selected) == count:
                return np.asarray(selected, dtype=int)
    raise ValueError(
        f"Found {len(selected)} variants at a minimum gap of {target_gap_bp} bp; needed {count}."
    )


def select_even_span(
    positions: np.ndarray,
    eligible_indices: np.ndarray,
    *,
    count: int,
    target_span_bp: int,
) -> np.ndarray:
    """Space `count` eligible variants evenly across `target_span_bp` from the first eligible position."""
    if count <= 0:
        raise ValueError("count must be positive.")
    if target_span_bp <= 0:
        raise ValueError("target_span_bp must be positive.")
    if len(eligible_indices) < count:
        raise ValueError(f"Eligible pool has {len(eligible_indices)} variants; needed {count}.")
    start_pos = int(positions[int(eligible_indices[0])])
    end_pos = start_pos + target_span_bp
    in_span = eligible_indices[positions[eligible_indices] <= end_pos]
    if len(in_span) < count:
        raise ValueError(
            f"Span of {target_span_bp} bp contains {len(in_span)} eligible variants; needed {count}."
        )
    choice = np.linspace(0, len(in_span) - 1, num=count)
    picked = np.unique(np.rint(choice).astype(int))
    if len(picked) != count:
        stride = len(in_span) / count
        picked = np.unique(np.minimum(len(in_span) - 1, np.floor(np.arange(count) * stride)).astype(int))
    if len(picked) != count:
        raise ValueError(f"Could not draw {count} unique variants from a {len(in_span)}-variant span.")
    return in_span[picked].astype(int, copy=False)


def greedy_ld_prune(
    dosage: np.ndarray,
    *,
    r2_threshold: float,
    count: int,
    chunk_size: int = 2048,
) -> np.ndarray:
    """Keep up to `count` variants in column order with pairwise r^2 at most `r2_threshold`.

    `dosage` is a samples-by-variants matrix already restricted to the candidate pool.
    Monomorphic columns are skipped. The result may be shorter than `count` when the
    pool cannot supply enough variants under the threshold.
    """
    if count <= 0:
        raise ValueError("count must be positive.")
    if not 0.0 <= r2_threshold <= 1.0:
        raise ValueError("r2_threshold must lie in [0, 1].")
    if dosage.ndim != 2:
        raise ValueError("dosage must be a samples-by-variants matrix.")
    n_samples, n_variants = dosage.shape
    if n_samples < 3 or n_variants == 0:
        raise ValueError("LD pruning needs at least 3 samples and 1 variant.")

    values = np.asarray(dosage, dtype=np.float64)
    mean = values.mean(axis=0)
    std = values.std(axis=0, ddof=1)
    variable = std > 1e-8
    standardized = np.zeros_like(values)
    standardized[:, variable] = (values[:, variable] - mean[variable]) / std[variable]
    denom = float(n_samples - 1)

    selected: list[int] = []
    selected_columns: list[np.ndarray] = []
    cursor = 0
    while cursor < n_variants and len(selected) < count:
        end = min(n_variants, cursor + chunk_size)
        block = standardized[:, cursor:end]
        local_count = end - cursor
        if selected_columns:
            selected_matrix = np.column_stack(selected_columns)
            correlation = (selected_matrix.T @ block) / denom
            max_r2 = np.max(correlation * correlation, axis=0)
        else:
            max_r2 = np.zeros(local_count, dtype=np.float64)
        gram = (block.T @ block) / denom
        accepted_local: list[int] = []
        for local in range(local_count):
            if len(selected) >= count:
                break
            column = cursor + local
            if not variable[column]:
                continue
            current = float(max_r2[local])
            if accepted_local:
                extra = gram[local, accepted_local]
                current = max(current, float(np.max(extra * extra)))
            if current <= r2_threshold + 1e-9:
                selected.append(column)
                accepted_local.append(local)
                selected_columns.append(standardized[:, column])
        cursor = end
    return np.asarray(selected, dtype=int)


def eligible_indices(
    positions: np.ndarray,
    *,
    start_index: int,
    end_pos: int | None,
    excluded: list[tuple[int, int]] | None = None,
) -> np.ndarray:
    """Variant indices at or after `start_index`, optionally inside a position cap and outside excluded intervals."""
    indices = np.arange(start_index, len(positions), dtype=int)
    if end_pos is not None:
        indices = indices[positions[indices] <= end_pos]
    if excluded and len(indices):
        keep = np.ones(len(indices), dtype=bool)
        pos = positions[indices]
        for start, end in excluded:
            keep &= ~((pos >= start) & (pos <= end))
        indices = indices[keep]
    return indices


def _overlaps(start_pos: int, end_pos: int, excluded: list[tuple[int, int]]) -> bool:
    return any(start_pos <= end and end_pos >= start for start, end in excluded)


def choose_anchor_indices(
    positions: np.ndarray,
    *,
    anchors: int,
    snp_count: int,
    min_gap_bp: int,
    excluded: list[tuple[int, int]],
    long_range_span_bp: int,
    sparse_gap_bp: int,
    ld_pool_span_bp: int = 30_000_000,
    min_following_variants: int | None = None,
) -> list[int]:
    """Pick start indices for additional regions that can host every regime.

    Anchors are spread across the chromosome, stay clear of excluded intervals
    for both the short dense block and the long-range span, and are separated
    by `min_gap_bp`. A chromosome-wide scan fills any band that has no valid start.
    """
    if anchors <= 0:
        raise ValueError("anchors must be positive.")
    n_variants = len(positions)
    if n_variants < snp_count:
        raise ValueError(f"Only {n_variants} variants; needed at least {snp_count}.")
    needed_following = snp_count if min_following_variants is None else min_following_variants

    def valid(start: int, chosen: list[int]) -> bool:
        if start < 0 or start + snp_count > n_variants:
            return False
        start_pos = int(positions[start])
        dense_end = int(positions[start + snp_count - 1])
        long_end = start_pos + long_range_span_bp
        sparse_end = start_pos + sparse_gap_bp * (snp_count - 1)
        if _overlaps(start_pos, dense_end, excluded):
            return False
        if _overlaps(start_pos, max(long_end, sparse_end), excluded):
            return False
        if long_end > int(positions[-1]) or sparse_end > int(positions[-1]):
            return False
        if any(abs(start_pos - int(positions[other])) < min_gap_bp for other in chosen):
            return False
        long_pool = eligible_indices(positions, start_index=start, end_pos=long_end, excluded=excluded)
        if len(long_pool) < snp_count:
            return False
        following = eligible_indices(
            positions,
            start_index=start,
            end_pos=start_pos + ld_pool_span_bp,
            excluded=excluded,
        )
        return len(following) >= needed_following

    selected: list[int] = []
    stride = max(snp_count, 200)

    def take_from(candidates: list[int]) -> None:
        for start in candidates:
            if len(selected) >= anchors:
                return
            if valid(start, selected):
                selected.append(start)
                return

    for band in range(anchors):
        lo = int(n_variants * band / anchors)
        hi = min(int(n_variants * (band + 1) / anchors), n_variants - snp_count + 1)
        if hi <= lo:
            continue
        mid = (lo + hi) // 2
        candidates = list(range(lo, hi, stride))
        candidates.sort(key=lambda index: (abs(index - mid), index))
        before = len(selected)
        take_from(candidates)
        if len(selected) == before:
            continue
    if len(selected) < anchors:
        fallback = list(range(0, n_variants - snp_count + 1, stride))
        for start in fallback:
            if len(selected) >= anchors:
                break
            if valid(start, selected):
                selected.append(start)
    if len(selected) != anchors:
        raise ValueError(f"Selected {len(selected)} anchors; requested {anchors}.")
    return selected


def panel_geometry(positions: np.ndarray) -> dict[str, float]:
    if len(positions) == 0:
        return {
            "start": 0,
            "end": 0,
            "span_bp": 0,
            "median_inter_snp_bp": float("nan"),
            "snps_per_mb": float("nan"),
        }
    ordered = np.sort(positions.astype(np.int64))
    start = int(ordered[0])
    end = int(ordered[-1])
    span = int(end - start)
    gaps = np.diff(ordered)
    median_gap = float(np.median(gaps)) if gaps.size else float("nan")
    snps_per_mb = float(len(ordered) / (span / 1_000_000)) if span > 0 else float("nan")
    return {
        "start": start,
        "end": end,
        "span_bp": span,
        "median_inter_snp_bp": median_gap,
        "snps_per_mb": snps_per_mb,
    }


def correlation_summaries(matrix: np.ndarray) -> dict[str, float]:
    """Mean absolute correlation and maximum pairwise r^2 across variant columns."""
    if matrix.shape[1] < 2:
        return {"mean_abs_corr": float("nan"), "max_pairwise_r2": float("nan")}
    values = matrix.astype(float)
    variable = values[:, np.std(values, axis=0) > 0.0]
    if variable.shape[1] < 2:
        return {"mean_abs_corr": 0.0, "max_pairwise_r2": 0.0}
    corr = np.corrcoef(variable, rowvar=False)
    upper = np.triu_indices(corr.shape[0], k=1)
    pair_corr = corr[upper]
    pair_corr = pair_corr[np.isfinite(pair_corr)]
    if pair_corr.size == 0:
        return {"mean_abs_corr": float("nan"), "max_pairwise_r2": float("nan")}
    r2 = pair_corr * pair_corr
    return {"mean_abs_corr": float(np.mean(np.abs(pair_corr))), "max_pairwise_r2": float(np.max(r2))}
