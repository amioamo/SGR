"""Reference haplotype reconstruction from aggregate Beacon carrier counts."""

from __future__ import annotations

from itertools import combinations_with_replacement
from typing import Any, Mapping

import numpy as np
from numpy.typing import ArrayLike

from .core import BeaconAttackProblem, ReconstructionResult
from .count_feasibility import assert_count_feasible

REFERENCE_HAPLOTYPE_RECONSTRUCTION = "reference_haplotype_reconstruction"
TEMPLATE_SELECTION_RULE = "proportional_greedy_with_replacement"
REFERENCE_ONLY = "reference_only"
RANDOM_TEMPLATE_PROJECTION = "random_template_projection"
GREEDY_TEMPLATE_WITHOUT_PROJECTION = "greedy_template_without_projection"
UNIFORM_ROW_SAMPLE = "uniform_row_sample"
EXACT_COUNT_PROJECTION = "exact_count_projection"


def _ensure_attack_problem(problem: Any, function_name: str) -> BeaconAttackProblem:
    if not isinstance(problem, BeaconAttackProblem):
        raise TypeError(f"{function_name} requires a BeaconAttackProblem.")
    return problem


def _as_binary_matrix(matrix: ArrayLike, name: str = "matrix") -> np.ndarray:
    array = np.asarray(matrix)
    if array.ndim != 2:
        raise ValueError(f"{name} must be two-dimensional.")
    if not np.all((array == 0) | (array == 1)):
        raise ValueError(f"{name} must contain only binary 0/1 values.")
    return array.astype(np.int8, copy=True)


def _resolve_seed(explicit_seed: int | None, configured_seed: int | None) -> int:
    if explicit_seed is not None:
        return int(explicit_seed)
    if configured_seed is not None:
        return int(configured_seed)
    return 0


def _json_safe(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _select_reference_templates(
    reference: np.ndarray,
    target_counts: np.ndarray,
    *,
    beacon_size: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Greedily track target frequencies with reusable reference rows."""

    if reference.shape[0] == 0:
        raise ValueError("reference_matrix must contain at least one row.")
    selected: list[int] = []
    current = np.zeros(reference.shape[1], dtype=np.int64)
    target = target_counts.astype(np.int64, copy=False)
    candidate_rows = reference.astype(np.int64, copy=False)
    for selected_count in range(1, beacon_size + 1):
        # Compare the partial template with the corresponding fraction of the
        # final target. Every row stays eligible, so row multiplicity acts as
        # a template weight even when the Beacon and reference sizes coincide.
        proposed_scaled_counts = beacon_size * (current[None, :] + candidate_rows)
        target_scaled_counts = selected_count * target[None, :]
        proposed_error = np.abs(target_scaled_counts - proposed_scaled_counts).sum(axis=1)
        best_error = int(np.min(proposed_error))
        best_positions = np.flatnonzero(proposed_error == best_error)
        chosen_index = int(rng.choice(best_positions))
        selected.append(chosen_index)
        current += candidate_rows[chosen_index]
    return reference[np.asarray(selected, dtype=int), :]


def minimum_template_count_l1(
    reference: ArrayLike,
    target_counts: ArrayLike,
    *,
    beacon_size: int,
    batch_size: int = 8192,
) -> int:
    """Minimum ``||y - c(T)||_1`` over multisets of ``beacon_size`` reference rows.

    Row reuse is allowed, matching Stage 1. Column sums do not depend on row
    order, so each multiset is enumerated once.
    """

    reference_rows = np.asarray(reference)
    if reference_rows.ndim != 2:
        raise ValueError("reference must be two-dimensional.")
    if reference_rows.shape[0] == 0:
        raise ValueError("reference must contain at least one row.")
    if beacon_size < 0:
        raise ValueError("beacon_size must be non-negative.")
    if batch_size < 1:
        raise ValueError("batch_size must be positive.")
    target = np.asarray(target_counts, dtype=np.int64)
    if target.shape != (reference_rows.shape[1],):
        raise ValueError("target_counts length must match reference columns.")
    if beacon_size == 0:
        return int(np.abs(target).sum())

    reference_rows = reference_rows.astype(np.int64, copy=False)
    row_count = int(reference_rows.shape[0])
    best: int | None = None
    batch: list[tuple[int, ...]] = []

    def fold(combos: list[tuple[int, ...]]) -> None:
        nonlocal best
        index = np.asarray(combos, dtype=np.intp)
        counts = np.zeros((index.shape[0], row_count), dtype=np.int64)
        np.add.at(counts, (np.arange(index.shape[0])[:, None], index), 1)
        column_sums = counts @ reference_rows
        l1 = np.abs(column_sums - target[None, :]).sum(axis=1)
        value = int(l1.min())
        best = value if best is None else min(best, value)

    for combo in combinations_with_replacement(range(row_count), beacon_size):
        batch.append(combo)
        if len(batch) >= batch_size:
            fold(batch)
            batch.clear()
    if batch:
        fold(batch)
    if best is None:
        raise RuntimeError("template enumeration produced no candidates.")
    return best


def _project_templates_to_counts(
    templates: np.ndarray,
    target_counts: np.ndarray,
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    beacon_size, variant_count = templates.shape
    matrix = np.zeros((beacon_size, variant_count), dtype=np.int8)
    random_tiebreaks = rng.random((beacon_size, variant_count))
    for column, count in enumerate(target_counts.astype(int)):
        if count <= 0:
            continue
        priority = templates[:, column].astype(float) + 1e-6 * random_tiebreaks[:, column]
        carrier_rows = np.argsort(priority)[-int(count) :]
        matrix[carrier_rows, column] = 1
    return matrix


def reference_haplotype_reconstruction(
    problem: BeaconAttackProblem,
    reference_matrix: ArrayLike,
    *,
    initialization_seed: int | None = None,
    optimizer_seed: int | None = None,
    config: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
) -> ReconstructionResult:
    """Reconstruct a count-feasible carrier matrix from reference haplotype rows."""

    attack_problem = _ensure_attack_problem(problem, "reference_haplotype_reconstruction")
    reference = _as_binary_matrix(reference_matrix, "reference_matrix")
    variant_count = len(attack_problem.variant_ids)
    if reference.shape[1] != variant_count:
        raise ValueError("reference_matrix column count must match problem variants.")

    target_counts = np.asarray(attack_problem.target_counts, dtype=np.int64)
    init_seed = _resolve_seed(initialization_seed, attack_problem.config.initialization_seed)
    opt_seed = _resolve_seed(optimizer_seed, attack_problem.config.optimizer_seed)
    rng = np.random.default_rng(init_seed)
    templates = _select_reference_templates(
        reference,
        target_counts,
        beacon_size=attack_problem.beacon_size,
        rng=rng,
    )
    reconstructed = _project_templates_to_counts(templates, target_counts, rng=rng)
    assert_count_feasible(target_counts, attack_problem.beacon_size, reconstructed)
    max_count_error = int(np.max(np.abs(reconstructed.sum(axis=0, dtype=np.int64) - target_counts), initial=0))
    seeds = {"initialization_seed": int(init_seed), "optimizer_seed": int(opt_seed)}
    if attack_problem.config.group_selection_seed is not None:
        seeds["group_selection_seed"] = int(attack_problem.config.group_selection_seed)
    default_config: dict[str, Any] = {
        "initialization_seed": int(init_seed),
        "optimizer_seed": int(opt_seed),
        "problem_config_metadata": _json_safe(dict(attack_problem.config.metadata)),
    }
    if config is not None:
        metadata_config = _json_safe(dict(config))
    elif metadata is not None:
        metadata_config = _json_safe(dict(metadata))
    else:
        metadata_config = default_config
    metadata_out: dict[str, Any] = {
        "method_name": REFERENCE_HAPLOTYPE_RECONSTRUCTION,
        "manifest_id": attack_problem.manifest_id,
        "split_id": attack_problem.split_id,
        "prior_id": "reference_haplotype_panel",
        "source_type": "reference_panel",
        "seeds": dict(seeds),
        "max_count_error": int(max_count_error),
        "initialization_source": "reference_row_selection_and_count_projection",
        "template_selection_rule": TEMPLATE_SELECTION_RULE,
        "config": metadata_config,
    }
    if metadata is not None and config is not None:
        metadata_out["run_metadata"] = _json_safe(dict(metadata))
    return ReconstructionResult(
        manifest_id=attack_problem.manifest_id,
        split_id=attack_problem.split_id,
        prior_id="reference_haplotype_panel",
        method_name=REFERENCE_HAPLOTYPE_RECONSTRUCTION,
        seeds=seeds,
        reconstructed_binary=reconstructed,
        metrics={"max_count_error": int(max_count_error)},
        metadata=metadata_out,
    )


def _baseline_seeds(attack_problem: BeaconAttackProblem, init_seed: int, opt_seed: int) -> dict[str, int]:
    seeds = {"initialization_seed": int(init_seed), "optimizer_seed": int(opt_seed)}
    if attack_problem.config.group_selection_seed is not None:
        seeds["group_selection_seed"] = int(attack_problem.config.group_selection_seed)
    return seeds


def _baseline_config(
    attack_problem: BeaconAttackProblem,
    init_seed: int,
    opt_seed: int,
    config: dict[str, Any] | None,
    metadata: dict[str, Any] | None,
) -> Any:
    default_config: dict[str, Any] = {
        "initialization_seed": int(init_seed),
        "optimizer_seed": int(opt_seed),
        "problem_config_metadata": _json_safe(dict(attack_problem.config.metadata)),
    }
    if config is not None:
        return _json_safe(dict(config))
    if metadata is not None:
        return _json_safe(dict(metadata))
    return default_config


def _prepare_reference_attack(
    problem: Any,
    reference_matrix: ArrayLike,
    function_name: str,
    initialization_seed: int | None,
    optimizer_seed: int | None,
) -> tuple[BeaconAttackProblem, np.ndarray, np.ndarray, dict[str, int], np.random.Generator, int, int]:
    attack_problem = _ensure_attack_problem(problem, function_name)
    reference = _as_binary_matrix(reference_matrix, "reference_matrix")
    if reference.shape[0] == 0:
        raise ValueError("reference_matrix must contain at least one row.")
    if reference.shape[1] != len(attack_problem.variant_ids):
        raise ValueError("reference_matrix column count must match problem variants.")
    target_counts = np.asarray(attack_problem.target_counts, dtype=np.int64)
    init_seed = _resolve_seed(initialization_seed, attack_problem.config.initialization_seed)
    opt_seed = _resolve_seed(optimizer_seed, attack_problem.config.optimizer_seed)
    seeds = _baseline_seeds(attack_problem, init_seed, opt_seed)
    rng = np.random.default_rng(init_seed)
    return attack_problem, reference, target_counts, seeds, rng, init_seed, opt_seed


def _sample_reference_rows(
    reference: np.ndarray,
    *,
    beacon_size: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, bool]:
    """Draw ``beacon_size`` reference rows without using released counts."""

    row_count = int(reference.shape[0])
    if row_count == 0:
        raise ValueError("reference_matrix must contain at least one row.")
    if beacon_size < 0:
        raise ValueError("beacon_size must be non-negative.")
    with_replacement = row_count < beacon_size
    indices = rng.choice(row_count, size=beacon_size, replace=with_replacement)
    return reference[np.asarray(indices, dtype=int), :], bool(with_replacement)


def _finish_stage_baseline(
    attack_problem: BeaconAttackProblem,
    reconstructed: np.ndarray,
    target_counts: np.ndarray,
    *,
    method_name: str,
    prior_id: str,
    seeds: dict[str, int],
    init_seed: int,
    opt_seed: int,
    initialization_source: str,
    template_selection_rule: str,
    stage1: str,
    stage2: str,
    uses_target_counts: bool,
    applies_count_projection: bool,
    sampling_with_replacement: bool | None,
    config: dict[str, Any] | None,
    metadata: dict[str, Any] | None,
) -> ReconstructionResult:
    max_count_error = int(np.max(np.abs(reconstructed.sum(axis=0, dtype=np.int64) - target_counts), initial=0))
    metadata_out: dict[str, Any] = {
        "method_name": method_name,
        "manifest_id": attack_problem.manifest_id,
        "split_id": attack_problem.split_id,
        "prior_id": prior_id,
        "source_type": "reference_panel",
        "seeds": dict(seeds),
        "max_count_error": int(max_count_error),
        "initialization_source": initialization_source,
        "template_selection_rule": template_selection_rule,
        "stage1": stage1,
        "stage2": stage2,
        "uses_target_counts": bool(uses_target_counts),
        "applies_count_projection": bool(applies_count_projection),
        "config": _baseline_config(attack_problem, init_seed, opt_seed, config, metadata),
    }
    if sampling_with_replacement is not None:
        metadata_out["sampling_with_replacement"] = bool(sampling_with_replacement)
    if metadata is not None and config is not None:
        metadata_out["run_metadata"] = _json_safe(dict(metadata))
    return ReconstructionResult(
        manifest_id=attack_problem.manifest_id,
        split_id=attack_problem.split_id,
        prior_id=prior_id,
        method_name=method_name,
        seeds=seeds,
        reconstructed_binary=reconstructed,
        metrics={"max_count_error": int(max_count_error)},
        metadata=metadata_out,
    )


def reference_only(
    problem: BeaconAttackProblem,
    reference_matrix: ArrayLike,
    *,
    initialization_seed: int | None = None,
    optimizer_seed: int | None = None,
    config: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
) -> ReconstructionResult:
    """Return N rows sampled from the reference panel, without using ``y``."""

    attack_problem, reference, target_counts, seeds, rng, init_seed, opt_seed = _prepare_reference_attack(
        problem,
        reference_matrix,
        "reference_only",
        initialization_seed,
        optimizer_seed,
    )
    reconstructed, with_replacement = _sample_reference_rows(
        reference,
        beacon_size=attack_problem.beacon_size,
        rng=rng,
    )
    return _finish_stage_baseline(
        attack_problem,
        reconstructed,
        target_counts,
        method_name=REFERENCE_ONLY,
        prior_id="reference_only_panel",
        seeds=seeds,
        init_seed=init_seed,
        opt_seed=opt_seed,
        initialization_source="uniform_reference_row_sample",
        template_selection_rule=UNIFORM_ROW_SAMPLE,
        stage1=UNIFORM_ROW_SAMPLE,
        stage2="none",
        uses_target_counts=False,
        applies_count_projection=False,
        sampling_with_replacement=with_replacement,
        config=config,
        metadata=metadata,
    )


def random_template_projection(
    problem: BeaconAttackProblem,
    reference_matrix: ArrayLike,
    *,
    initialization_seed: int | None = None,
    optimizer_seed: int | None = None,
    config: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
) -> ReconstructionResult:
    """Sample N reference rows, then apply the Stage 2 exact count projection."""

    attack_problem, reference, target_counts, seeds, rng, init_seed, opt_seed = _prepare_reference_attack(
        problem,
        reference_matrix,
        "random_template_projection",
        initialization_seed,
        optimizer_seed,
    )
    templates, with_replacement = _sample_reference_rows(
        reference,
        beacon_size=attack_problem.beacon_size,
        rng=rng,
    )
    reconstructed = _project_templates_to_counts(templates, target_counts, rng=rng)
    assert_count_feasible(target_counts, attack_problem.beacon_size, reconstructed)
    return _finish_stage_baseline(
        attack_problem,
        reconstructed,
        target_counts,
        method_name=RANDOM_TEMPLATE_PROJECTION,
        prior_id="random_template_panel",
        seeds=seeds,
        init_seed=init_seed,
        opt_seed=opt_seed,
        initialization_source="uniform_reference_row_sample_and_count_projection",
        template_selection_rule=UNIFORM_ROW_SAMPLE,
        stage1=UNIFORM_ROW_SAMPLE,
        stage2=EXACT_COUNT_PROJECTION,
        uses_target_counts=True,
        applies_count_projection=True,
        sampling_with_replacement=with_replacement,
        config=config,
        metadata=metadata,
    )


def greedy_template_without_projection(
    problem: BeaconAttackProblem,
    reference_matrix: ArrayLike,
    *,
    initialization_seed: int | None = None,
    optimizer_seed: int | None = None,
    config: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
) -> ReconstructionResult:
    """Return the Stage 1 greedy template without the Stage 2 count projection."""

    attack_problem, reference, target_counts, seeds, rng, init_seed, opt_seed = _prepare_reference_attack(
        problem,
        reference_matrix,
        "greedy_template_without_projection",
        initialization_seed,
        optimizer_seed,
    )
    reconstructed = _select_reference_templates(
        reference,
        target_counts,
        beacon_size=attack_problem.beacon_size,
        rng=rng,
    )
    return _finish_stage_baseline(
        attack_problem,
        reconstructed,
        target_counts,
        method_name=GREEDY_TEMPLATE_WITHOUT_PROJECTION,
        prior_id="greedy_template_panel",
        seeds=seeds,
        init_seed=init_seed,
        opt_seed=opt_seed,
        initialization_source="greedy_reference_template",
        template_selection_rule=TEMPLATE_SELECTION_RULE,
        stage1=TEMPLATE_SELECTION_RULE,
        stage2="none",
        uses_target_counts=True,
        applies_count_projection=False,
        sampling_with_replacement=None,
        config=config,
        metadata=metadata,
    )
