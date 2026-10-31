"""Shared helpers for experiment timing and convergence metadata."""

from __future__ import annotations

import time
from collections.abc import Callable
from typing import Any, TypeVar

from beacon_reconstruction.core import ReconstructionResult

T = TypeVar("T")


def run_timed(fn: Callable[[], T]) -> tuple[T, float]:
    start = time.perf_counter()
    result = fn()
    return result, time.perf_counter() - start


def runtime_metadata_fields(result: ReconstructionResult, wall_time_sec: float) -> dict[str, Any]:
    fields: dict[str, Any] = {
        "wall_time_sec": round(float(wall_time_sec), 6),
        "reference_panel_mode": result.metadata.get("reference_panel_mode", ""),
    }
    meta = result.metadata or {}
    if "outer_iterations" in meta:
        fields["saleem_outer_iterations"] = int(meta["outer_iterations"])
        fields["saleem_converged"] = bool(meta.get("converged", False))
    if "flip_count" in meta:
        fields["saleem_flip_count"] = int(meta["flip_count"])
    loss_trace = meta.get("loss_trace")
    if isinstance(loss_trace, list) and loss_trace:
        fields["pairwise_epochs_executed"] = max(0, len(loss_trace) - 1)
    if "accepted_moves" in meta:
        fields["pairwise_accepted_moves"] = int(meta["accepted_moves"])
    config = meta.get("config")
    if isinstance(config, dict):
        nested = config.get("config")
        if isinstance(nested, dict) and nested.get("pairwise_epochs") is not None:
            fields["pairwise_epochs_budget"] = int(nested["pairwise_epochs"])
        elif config.get("pairwise_epochs") is not None:
            fields["pairwise_epochs_budget"] = int(config["pairwise_epochs"])
    return fields


def attach_reference_panel_mode(result: ReconstructionResult, *, shared: bool) -> ReconstructionResult:
    metadata = dict(result.metadata)
    metadata["reference_panel_mode"] = "shared" if shared else "dual_disjoint"
    return ReconstructionResult(
        manifest_id=result.manifest_id,
        split_id=result.split_id,
        prior_id=result.prior_id,
        method_name=result.method_name,
        seeds=result.seeds,
        reconstructed_binary=result.reconstructed_binary,
        metrics=result.metrics,
        metadata=metadata,
    )
