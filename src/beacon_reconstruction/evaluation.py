"""Evaluation utilities for Beacon reconstruction experiments."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import permutations
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .core import BeaconEvaluationTruth, ReconstructionResult

try:  # pragma: no cover
    from scipy.optimize import linear_sum_assignment
except ImportError:  # pragma: no cover
    linear_sum_assignment = None


ALLOWED_METHOD_NAMES = frozenset(
    {
        "independent_count",
        "pairwise_refinement",
        "saleem_reconstruction_attack",
        "reference_haplotype_reconstruction",
        "reference_only",
        "random_template_projection",
        "greedy_template_without_projection",
    }
)


@dataclass(frozen=True)
class RowAlignment:
    """One-to-one assignment between true and reconstructed rows."""

    true_indices: NDArray[np.integer[Any]]
    reconstructed_indices: NDArray[np.integer[Any]]
    total_hamming_cost: int
    cost_matrix: NDArray[np.integer[Any]]


@dataclass(frozen=True)
class ClassificationMetrics:
    """Binary reconstruction metrics after row alignment."""

    precision: float
    recall: float
    f1: float
    true_positives: int
    false_positives: int
    false_negatives: int
    alignment_total_hamming_cost: int


@dataclass(frozen=True)
class CountMetrics:
    """Carrier-count feasibility metrics."""

    count_error: NDArray[np.integer[Any]]
    count_mae: float
    max_count_error: float


@dataclass(frozen=True)
class EvaluationMetrics:
    """Per-run raw metrics. NaN values are preserved."""

    precision: float
    recall: float
    f1: float
    count_error: NDArray[np.integer[Any]]
    count_mae: float
    max_count_error: float
    pairwise_error: float
    true_positives: int
    false_positives: int
    false_negatives: int
    alignment_total_hamming_cost: int

    def as_dict(self) -> dict[str, Any]:
        return {
            "precision": self.precision,
            "recall": self.recall,
            "f1": self.f1,
            "count_error": self.count_error.copy(),
            "count_mae": self.count_mae,
            "max_count_error": self.max_count_error,
            "pairwise_error": self.pairwise_error,
            "true_positives": self.true_positives,
            "false_positives": self.false_positives,
            "false_negatives": self.false_negatives,
            "alignment_total_hamming_cost": self.alignment_total_hamming_cost,
        }


@dataclass(frozen=True)
class InvalidRunReport:
    """Invalid-run status plus recorded reasons."""

    is_invalid: bool
    reasons: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class EvaluationResult:
    """Evaluation metrics bundled with invalid-run status."""

    metrics: EvaluationMetrics
    invalid_run: InvalidRunReport


def _as_binary_matrix(values: ArrayLike, name: str) -> NDArray[np.int8]:
    matrix = np.asarray(values)
    if matrix.ndim != 2:
        raise ValueError(f"{name} must be a two-dimensional matrix.")
    if not np.all((matrix == 0) | (matrix == 1)):
        raise ValueError(f"{name} must be binary.")
    return matrix.astype(np.int8, copy=True)


def _variant_subset(matrix: NDArray[np.int8], evaluated_variant_indices: Sequence[int] | None) -> NDArray[np.int8]:
    if evaluated_variant_indices is None:
        return matrix
    indices = np.asarray(evaluated_variant_indices, dtype=int)
    if indices.ndim != 1:
        raise ValueError("evaluated_variant_indices must be one-dimensional.")
    if np.any(indices < 0) or np.any(indices >= matrix.shape[1]):
        raise ValueError("evaluated_variant_indices contains an out-of-bounds variant index.")
    return matrix[:, indices]


def _safe_divide(numerator: int | float, denominator: int | float) -> float:
    if denominator == 0:
        return float("nan")
    return float(numerator / denominator)


def _hungarian_assignment(cost_matrix: NDArray[np.integer[Any]]) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    if linear_sum_assignment is not None:
        row_indices, column_indices = linear_sum_assignment(cost_matrix)
        return row_indices.astype(np.int64), column_indices.astype(np.int64)

    row_count, column_count = cost_matrix.shape
    if row_count != column_count:
        raise RuntimeError("scipy is required for rectangular Hungarian matching.")
    best_columns: tuple[int, ...] | None = None
    best_cost: int | None = None
    for columns in permutations(range(column_count)):
        cost = int(sum(cost_matrix[row, column] for row, column in enumerate(columns)))
        if best_cost is None or cost < best_cost:
            best_cost = cost
            best_columns = columns
    assert best_columns is not None and best_cost is not None
    return np.arange(row_count, dtype=np.int64), np.asarray(best_columns, dtype=np.int64)


def align_rows_hungarian(
    true_binary: ArrayLike,
    reconstructed_binary: ArrayLike,
    *,
    evaluated_variant_indices: Sequence[int] | None = None,
) -> RowAlignment:
    """Align unordered reconstructed rows to truth using one-to-one Hungarian matching."""

    true_matrix = _as_binary_matrix(true_binary, "true_binary")
    reconstructed_matrix = _as_binary_matrix(reconstructed_binary, "reconstructed_binary")
    if true_matrix.shape != reconstructed_matrix.shape:
        raise ValueError("true_binary and reconstructed_binary must have the same shape.")

    true_eval = _variant_subset(true_matrix, evaluated_variant_indices)
    reconstructed_eval = _variant_subset(reconstructed_matrix, evaluated_variant_indices)
    cost_matrix = np.not_equal(true_eval[:, None, :], reconstructed_eval[None, :, :]).sum(axis=2)
    true_indices, reconstructed_indices = _hungarian_assignment(cost_matrix)
    if len(np.unique(true_indices)) != len(true_indices) or len(np.unique(reconstructed_indices)) != len(
        reconstructed_indices
    ):
        raise RuntimeError("Hungarian matching failed to produce a one-to-one assignment.")
    total_hamming_cost = int(cost_matrix[true_indices, reconstructed_indices].sum())
    return RowAlignment(
        true_indices=true_indices,
        reconstructed_indices=reconstructed_indices,
        total_hamming_cost=total_hamming_cost,
        cost_matrix=cost_matrix.astype(int, copy=True),
    )


def compute_classification_metrics_after_alignment(
    true_binary: ArrayLike,
    reconstructed_binary: ArrayLike,
    *,
    evaluated_variant_indices: Sequence[int] | None = None,
) -> ClassificationMetrics:
    """Compute precision, recall, and F1 after one-to-one row alignment."""

    true_matrix = _as_binary_matrix(true_binary, "true_binary")
    reconstructed_matrix = _as_binary_matrix(reconstructed_binary, "reconstructed_binary")
    alignment = align_rows_hungarian(
        true_matrix,
        reconstructed_matrix,
        evaluated_variant_indices=evaluated_variant_indices,
    )
    true_eval = _variant_subset(true_matrix, evaluated_variant_indices)[alignment.true_indices]
    reconstructed_eval = _variant_subset(reconstructed_matrix, evaluated_variant_indices)[alignment.reconstructed_indices]

    true_positives = int(np.logical_and(true_eval == 1, reconstructed_eval == 1).sum())
    false_positives = int(np.logical_and(true_eval == 0, reconstructed_eval == 1).sum())
    false_negatives = int(np.logical_and(true_eval == 1, reconstructed_eval == 0).sum())

    precision = _safe_divide(true_positives, true_positives + false_positives)
    recall = _safe_divide(true_positives, true_positives + false_negatives)
    f1 = float("nan")
    if np.isfinite(precision) and np.isfinite(recall) and (precision + recall) != 0:
        f1 = float(2 * precision * recall / (precision + recall))

    return ClassificationMetrics(
        precision=precision,
        recall=recall,
        f1=f1,
        true_positives=true_positives,
        false_positives=false_positives,
        false_negatives=false_negatives,
        alignment_total_hamming_cost=alignment.total_hamming_cost,
    )


def compute_count_metrics(reconstructed_binary: ArrayLike, target_counts: ArrayLike) -> CountMetrics:
    """Compute count errors against carrier-count targets."""

    reconstructed_matrix = _as_binary_matrix(reconstructed_binary, "reconstructed_binary")
    counts = np.asarray(target_counts)
    if counts.ndim != 1:
        raise ValueError("target_counts must be one-dimensional.")
    if len(counts) != reconstructed_matrix.shape[1]:
        raise ValueError("target_counts length must match reconstructed_binary columns.")
    if not np.issubdtype(counts.dtype, np.integer):
        if not np.all(np.equal(counts, np.floor(counts))):
            raise ValueError("target_counts must be integer carrier counts.")
        counts = counts.astype(int)
    count_error = reconstructed_matrix.sum(axis=0).astype(int) - counts.astype(int)
    if count_error.size == 0:
        count_mae = float("nan")
        max_count_error = float("nan")
    else:
        count_mae = float(np.mean(np.abs(count_error)))
        max_count_error = float(np.max(np.abs(count_error)))
    return CountMetrics(count_error=count_error.astype(int, copy=True), count_mae=count_mae, max_count_error=max_count_error)


def compute_pairwise_matrix(binary_matrix: ArrayLike) -> NDArray[np.float64]:
    """Compute empirical pairwise co-carrier probabilities."""

    matrix = _as_binary_matrix(binary_matrix, "binary_matrix").astype(float)
    if matrix.shape[0] == 0:
        return np.full((matrix.shape[1], matrix.shape[1]), np.nan, dtype=float)
    return (matrix.T @ matrix) / matrix.shape[0]


def compute_pairwise_error(
    reconstructed_binary: ArrayLike,
    reference_pairwise: ArrayLike | None,
    *,
    evaluated_variant_indices: Sequence[int] | None = None,
) -> float:
    """Mean squared error over upper-triangle pairwise co-carrier moments."""

    if reference_pairwise is None:
        return float("nan")
    reconstructed_matrix = _as_binary_matrix(reconstructed_binary, "reconstructed_binary")
    reconstructed_eval = _variant_subset(reconstructed_matrix, evaluated_variant_indices)
    observed = compute_pairwise_matrix(reconstructed_eval)
    reference = np.asarray(reference_pairwise, dtype=float)
    if evaluated_variant_indices is not None and reference.shape == (
        reconstructed_matrix.shape[1],
        reconstructed_matrix.shape[1],
    ):
        indices = np.asarray(evaluated_variant_indices, dtype=int)
        reference = reference[np.ix_(indices, indices)]
    if reference.shape != observed.shape:
        raise ValueError("reference_pairwise must match the evaluated variant count.")
    if observed.shape[0] < 2:
        return float("nan")
    upper = np.triu_indices(observed.shape[0], k=1)
    return float(np.mean((observed[upper] - reference[upper]) ** 2))


def evaluate_reconstruction(
    result: ReconstructionResult,
    truth: BeaconEvaluationTruth,
    target_counts: ArrayLike,
    *,
    reference_pairwise: ArrayLike | None = None,
    evaluated_variant_indices: Sequence[int] | None = None,
    source_sets: Mapping[str, Sequence[str]] | None = None,
    accepted_manifest_ids: Sequence[str] | None = None,
) -> EvaluationResult:
    """Evaluate a completed reconstruction and record invalid-run status."""

    if result.manifest_id != truth.manifest_id or result.split_id != truth.split_id:
        raise ValueError("result and truth must have matching manifest_id and split_id.")

    reconstructed = result.reconstructed_binary
    true_binary = truth.true_beacon_binary
    if evaluated_variant_indices is not None:
        target_counts_array = np.asarray(target_counts)
        target_counts = target_counts_array[np.asarray(evaluated_variant_indices, dtype=int)]
        reconstructed_for_counts = _variant_subset(reconstructed, evaluated_variant_indices)
    else:
        reconstructed_for_counts = reconstructed

    classification = compute_classification_metrics_after_alignment(
        true_binary,
        reconstructed,
        evaluated_variant_indices=evaluated_variant_indices,
    )
    count_metrics = compute_count_metrics(reconstructed_for_counts, target_counts)
    pairwise_error = compute_pairwise_error(
        reconstructed,
        reference_pairwise,
        evaluated_variant_indices=evaluated_variant_indices,
    )

    metrics = EvaluationMetrics(
        precision=classification.precision,
        recall=classification.recall,
        f1=classification.f1,
        count_error=count_metrics.count_error,
        count_mae=count_metrics.count_mae,
        max_count_error=count_metrics.max_count_error,
        pairwise_error=pairwise_error,
        true_positives=classification.true_positives,
        false_positives=classification.false_positives,
        false_negatives=classification.false_negatives,
        alignment_total_hamming_cost=classification.alignment_total_hamming_cost,
    )
    invalid_run = check_invalid_run(
        result,
        metrics=metrics,
        source_sets=source_sets,
        accepted_manifest_ids=accepted_manifest_ids,
    )
    return EvaluationResult(metrics=metrics, invalid_run=invalid_run)


def _get_field(record: Any, name: str) -> Any:
    if isinstance(record, Mapping):
        return record.get(name)
    return getattr(record, name, None)


def _metadata_from(record: Any) -> Mapping[str, Any]:
    metadata = _get_field(record, "metadata")
    if isinstance(metadata, Mapping):
        return metadata
    return {}


def _is_present(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str) and not value:
        return False
    if isinstance(value, Mapping) and len(value) == 0:
        return False
    return True


def _used_evaluation_truth(metadata: Mapping[str, Any]) -> bool:
    flags = (
        "used_evaluation_truth",
        "method_used_evaluation_truth",
        "evaluation_truth_used",
        "truth_leakage",
    )
    return any(bool(metadata.get(flag, False)) for flag in flags)


def _source_overlap_reasons(source_sets: Mapping[str, Sequence[str]] | None) -> list[str]:
    if source_sets is None:
        return []
    normalized = {name: set(values) for name, values in source_sets.items()}
    beacon_names = [name for name in normalized if "beacon" in name.lower()]
    reference_names = [
        name
        for name in normalized
        if name not in beacon_names and ("reference" in name.lower() or "ref" in name.lower() or "arg" in name.lower())
    ]
    reasons: list[str] = []
    for beacon_name in beacon_names:
        for reference_name in reference_names:
            overlap = normalized[beacon_name] & normalized[reference_name]
            if overlap:
                reasons.append(f"source_overlap:{beacon_name}:{reference_name}:{len(overlap)}")
    return reasons


def check_invalid_run(
    result: Any,
    *,
    metrics: EvaluationMetrics | Mapping[str, Any] | None = None,
    source_sets: Mapping[str, Sequence[str]] | None = None,
    allowed_method_names: frozenset[str] = ALLOWED_METHOD_NAMES,
    accepted_manifest_ids: Sequence[str] | None = None,
) -> InvalidRunReport:
    """Check invalid-run rules without deleting or mutating the run."""

    reasons: list[str] = []
    for field_name in ("manifest_id", "split_id", "prior_id", "method_name"):
        if not _is_present(_get_field(result, field_name)):
            reasons.append(f"missing_{field_name}")

    manifest_id = _get_field(result, "manifest_id")
    if (
        accepted_manifest_ids is not None
        and _is_present(manifest_id)
        and str(manifest_id) not in set(accepted_manifest_ids)
    ):
        reasons.append("manifest_id_not_accepted")

    method_name = _get_field(result, "method_name")
    if _is_present(method_name) and method_name not in allowed_method_names:
        reasons.append("method_name_not_allowed")

    seeds = _get_field(result, "seeds")
    if not isinstance(seeds, Mapping) or not seeds:
        reasons.append("missing_seeds")
    else:
        for seed_name, seed_value in seeds.items():
            if not isinstance(seed_name, str) or not seed_name:
                reasons.append("invalid_seed_name")
            elif not isinstance(seed_value, int):
                reasons.append(f"invalid_seed:{seed_name}")

    metadata = _metadata_from(result)
    if not _is_present(metadata.get("config")):
        reasons.append("missing_config")

    if _used_evaluation_truth(metadata):
        reasons.append("used_evaluation_truth")
    if bool(metadata.get("source_overlap", False)):
        reasons.append("source_overlap")
    reasons.extend(_source_overlap_reasons(source_sets))

    max_count_error: Any = None
    if metrics is not None:
        max_count_error = metrics["max_count_error"] if isinstance(metrics, Mapping) else metrics.max_count_error
    else:
        result_metrics = _get_field(result, "metrics")
        if isinstance(result_metrics, Mapping):
            max_count_error = result_metrics.get("max_count_error")
    if max_count_error is not None:
        try:
            max_count_error_value = float(max_count_error)
        except (TypeError, ValueError):
            reasons.append("invalid_max_count_error")
        else:
            if np.isfinite(max_count_error_value) and max_count_error_value > 0:
                reasons.append("count_violation")

    return InvalidRunReport(is_invalid=bool(reasons), reasons=tuple(reasons))


@dataclass(frozen=True)
class PartialGenotypeLinkageMetrics:
    """Linkage accuracy from partial-genotype nearest-neighbor matching."""

    known_snp_count: int
    n_participants: int
    top1_accuracy: float
    top5_accuracy: float
    mean_rank_of_true_match: float
    median_rank_of_true_match: float


def select_partial_genotype_indices(
    variant_count: int,
    known_snp_count: int,
    *,
    seed: int,
) -> NDArray[np.integer[Any]]:
    """Choose a deterministic subset of variant indices for a partial-genotype query."""

    if known_snp_count < 0:
        raise ValueError("known_snp_count must be non-negative.")
    if known_snp_count > variant_count:
        raise ValueError("known_snp_count cannot exceed variant_count.")
    if known_snp_count == 0:
        return np.asarray([], dtype=int)
    rng = np.random.default_rng(int(seed))
    return np.sort(rng.choice(variant_count, size=known_snp_count, replace=False).astype(int))


def compute_partial_genotype_linkage(
    true_binary: ArrayLike,
    reconstructed_binary: ArrayLike,
    known_snp_count: int,
    *,
    base_seed: int,
    alignment: RowAlignment | None = None,
    top_k: int = 5,
) -> PartialGenotypeLinkageMetrics:
    """Match partial genotypes to reconstructed rows and report top-k linkage accuracy.

    After anonymous one-to-one row alignment, each true Beacon participant is paired
    with a reconstructed row. For each participant, an attacker who knows carrier
    status at ``known_snp_count`` variants ranks reconstructed rows by Hamming
    distance on those variants and checks whether the aligned row is near the top.
    """

    true_matrix = _as_binary_matrix(true_binary, "true_binary")
    reconstructed_matrix = _as_binary_matrix(reconstructed_binary, "reconstructed_binary")
    if true_matrix.shape != reconstructed_matrix.shape:
        raise ValueError("true_binary and reconstructed_binary must have the same shape.")
    if known_snp_count < 0:
        raise ValueError("known_snp_count must be non-negative.")
    if top_k < 1:
        raise ValueError("top_k must be at least 1.")

    resolved_alignment = alignment or align_rows_hungarian(true_matrix, reconstructed_matrix)
    variant_count = int(true_matrix.shape[1])
    participant_count = int(true_matrix.shape[0])
    effective_top_k = min(int(top_k), participant_count)

    top1_hits = 0
    top5_hits = 0
    ranks: list[int] = []
    for true_index, reconstructed_index in zip(
        resolved_alignment.true_indices,
        resolved_alignment.reconstructed_indices,
        strict=True,
    ):
        known_indices = select_partial_genotype_indices(
            variant_count,
            known_snp_count,
            seed=int(base_seed) + int(true_index) * 1009 + known_snp_count,
        )
        query = true_matrix[int(true_index), known_indices]
        candidate_distances = np.not_equal(
            query[None, :],
            reconstructed_matrix[:, known_indices],
        ).sum(axis=1)
        ranked_indices = np.argsort(candidate_distances, kind="stable")
        rank = int(np.where(ranked_indices == int(reconstructed_index))[0][0]) + 1
        ranks.append(rank)
        if rank == 1:
            top1_hits += 1
        if rank <= effective_top_k:
            top5_hits += 1

    rank_array = np.asarray(ranks, dtype=float)
    return PartialGenotypeLinkageMetrics(
        known_snp_count=known_snp_count,
        n_participants=participant_count,
        top1_accuracy=float(top1_hits / participant_count) if participant_count else float("nan"),
        top5_accuracy=float(top5_hits / participant_count) if participant_count else float("nan"),
        mean_rank_of_true_match=float(np.mean(rank_array)) if rank_array.size else float("nan"),
        median_rank_of_true_match=float(np.median(rank_array)) if rank_array.size else float("nan"),
    )
