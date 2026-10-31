"""Core objects for Beacon carrier-count reconstruction."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray

ResponseChannel = Literal["carrier_count"]


def _as_1d_array(values: ArrayLike, name: str, dtype: Any | None = None) -> NDArray[Any]:
    array = np.asarray(values, dtype=dtype)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array.")
    return array


def _as_optional_2d_float(values: ArrayLike | None, name: str) -> NDArray[np.float64] | None:
    if values is None:
        return None
    array = np.asarray(values, dtype=float)
    if array.ndim != 2:
        raise ValueError(f"{name} must be a two-dimensional array.")
    if np.any(~np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values.")
    return array


def _as_string_tuple(values: ArrayLike, name: str) -> tuple[str, ...]:
    array = _as_1d_array(values, name)
    return tuple(str(value) for value in array.tolist())


@dataclass(frozen=True)
class ExperimentConfig:
    """Run configuration passed to reconstruction methods."""

    response_channel: ResponseChannel = "carrier_count"
    initialization_seed: int | None = None
    optimizer_seed: int | None = None
    group_selection_seed: int | None = None
    random_seed: int | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.response_channel != "carrier_count":
            raise ValueError("Only response_channel='carrier_count' is supported.")


@dataclass(frozen=True)
class BeaconAttackProblem:
    """Attack-time inputs for carrier-count reconstruction."""

    manifest_id: str
    split_id: str
    beacon_size: int
    variant_ids: tuple[str, ...]
    positions: tuple[int, ...]
    query_responses: NDArray[np.integer[Any]]
    response_channel: ResponseChannel = "carrier_count"
    target_counts: NDArray[np.integer[Any]] | None = None
    public_ld: NDArray[np.float64] | None = None
    config: ExperimentConfig = field(default_factory=ExperimentConfig)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.manifest_id:
            raise ValueError("manifest_id is required.")
        if not self.split_id:
            raise ValueError("split_id is required.")
        if self.beacon_size < 0:
            raise ValueError("beacon_size must be non-negative.")
        if self.response_channel != "carrier_count":
            raise ValueError("Only response_channel='carrier_count' is supported.")

        variant_ids = _as_string_tuple(self.variant_ids, "variant_ids")
        positions_array = _as_1d_array(self.positions, "positions")
        positions = tuple(int(position) for position in positions_array.tolist())
        query_responses = _as_1d_array(self.query_responses, "query_responses")
        if len(variant_ids) != len(positions) or len(variant_ids) != len(query_responses):
            raise ValueError("variant_ids, positions, and query_responses must have the same length.")
        if not np.issubdtype(query_responses.dtype, np.integer):
            if not np.all(np.equal(query_responses, np.floor(query_responses))):
                raise ValueError("query_responses must be integer carrier counts.")
            query_responses = query_responses.astype(int)

        if self.target_counts is None:
            target_counts = query_responses.astype(int, copy=True)
        else:
            target_counts = _as_1d_array(self.target_counts, "target_counts")
            if len(target_counts) != len(query_responses):
                raise ValueError("target_counts must match query_responses length.")
            if not np.issubdtype(target_counts.dtype, np.integer):
                if not np.all(np.equal(target_counts, np.floor(target_counts))):
                    raise ValueError("target_counts must be integer carrier counts.")
                target_counts = target_counts.astype(int)
            if not np.array_equal(target_counts, query_responses):
                raise ValueError("For carrier_count responses, target_counts must equal query_responses.")

        if np.any(target_counts < 0) or np.any(target_counts > self.beacon_size):
            raise ValueError("target_counts must be between 0 and beacon_size.")

        public_ld = _as_optional_2d_float(self.public_ld, "public_ld")
        if public_ld is not None and public_ld.shape != (len(variant_ids), len(variant_ids)):
            raise ValueError("public_ld must have shape (variant_count, variant_count).")

        object.__setattr__(self, "variant_ids", variant_ids)
        object.__setattr__(self, "positions", positions)
        object.__setattr__(self, "query_responses", query_responses.astype(int, copy=True))
        object.__setattr__(self, "target_counts", target_counts.astype(int, copy=True))
        object.__setattr__(self, "public_ld", public_ld)


@dataclass(frozen=True)
class BeaconEvaluationTruth:
    """Evaluation-only fields that reconstruction functions must not consume."""

    manifest_id: str
    split_id: str
    true_beacon_binary: NDArray[np.integer[Any]]
    true_sample_ids: tuple[str, ...] = field(default_factory=tuple)
    true_row_order: tuple[int, ...] | None = None
    true_ancestry_labels: tuple[str, ...] = field(default_factory=tuple)
    evaluation_labels: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        truth = np.asarray(self.true_beacon_binary)
        if truth.ndim != 2:
            raise ValueError("true_beacon_binary must be a two-dimensional matrix.")
        if not np.all((truth == 0) | (truth == 1)):
            raise ValueError("true_beacon_binary must be binary.")
        object.__setattr__(self, "true_beacon_binary", truth.astype(np.int8, copy=True))
        object.__setattr__(self, "true_sample_ids", _as_string_tuple(self.true_sample_ids, "true_sample_ids"))
        object.__setattr__(
            self,
            "true_ancestry_labels",
            _as_string_tuple(self.true_ancestry_labels, "true_ancestry_labels"),
        )
        if self.true_row_order is not None:
            object.__setattr__(self, "true_row_order", tuple(int(i) for i in self.true_row_order))


@dataclass(frozen=True)
class ReconstructionResult:
    """Schema for a completed reconstruction output."""

    manifest_id: str
    split_id: str
    prior_id: str
    method_name: str
    seeds: Mapping[str, int]
    reconstructed_binary: NDArray[np.integer[Any]]
    metrics: Mapping[str, Any] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for field_name in ("manifest_id", "split_id", "prior_id", "method_name"):
            if not getattr(self, field_name):
                raise ValueError(f"{field_name} is required.")
        if not isinstance(self.seeds, Mapping) or not self.seeds:
            raise ValueError("seeds must record all relevant run seeds.")
        for seed_name, seed_value in self.seeds.items():
            if not isinstance(seed_name, str) or not seed_name:
                raise ValueError("seed names must be non-empty strings.")
            if not isinstance(seed_value, int):
                raise ValueError("seed values must be integers.")
        matrix = np.asarray(self.reconstructed_binary)
        if matrix.ndim != 2:
            raise ValueError("reconstructed_binary must be a two-dimensional matrix.")
        if not np.all((matrix == 0) | (matrix == 1)):
            raise ValueError("reconstructed_binary must be binary.")
        object.__setattr__(self, "reconstructed_binary", matrix.astype(np.int8, copy=True))
