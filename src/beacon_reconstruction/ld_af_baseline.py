"""LD/AF optimization baseline (Saleem et al., 2025).

Adapted from the official Beacon reconstruction attack implementation
(https://github.com/ASAP-Bilkent/Beacon-Reconstruction-Attack).
Matrix layout follows the original code: SNPs × individuals.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .core import BeaconAttackProblem, ReconstructionResult

try:
    import torch
except ImportError:  # pragma: no cover
    torch = None  # type: ignore[assignment]

SALEEM_RECONSTRUCTION_ATTACK = "saleem_reconstruction_attack"

__all__ = ["SALEEM_RECONSTRUCTION_ATTACK", "saleem_reconstruction_attack"]


def _ensure_attack_problem(problem: Any, function_name: str) -> BeaconAttackProblem:
    if not isinstance(problem, BeaconAttackProblem):
        raise TypeError(f"{function_name} requires a BeaconAttackProblem.")
    return problem


def _as_binary_matrix(binary_matrix: ArrayLike, name: str = "binary_matrix") -> NDArray[np.int8]:
    matrix = np.asarray(binary_matrix)
    if matrix.ndim != 2:
        raise ValueError(f"{name} must be a two-dimensional matrix.")
    if not np.all((matrix == 0) | (matrix == 1)):
        raise ValueError(f"{name} must contain only 0/1 values.")
    return matrix.astype(np.int8, copy=True)


def _resolve_seed(explicit_seed: int | None, configured_seed: int | None) -> int:
    if explicit_seed is not None:
        return int(explicit_seed)
    if configured_seed is not None:
        return int(configured_seed)
    return 0


def _require_torch() -> Any:
    if torch is None:
        raise ImportError(
            "saleem_reconstruction_attack requires PyTorch. Install with: pip install torch"
        )
    return torch


def _to_saleem_orientation(participant_matrix: ArrayLike) -> NDArray[np.int8]:
    matrix = _as_binary_matrix(participant_matrix)
    return matrix.T.copy()


def _from_saleem_orientation(saleem_matrix: ArrayLike) -> NDArray[np.int8]:
    matrix = _as_binary_matrix(saleem_matrix)
    return matrix.T.copy()


def _reference_panel_to_sokal_michener(reference_matrix: ArrayLike) -> NDArray[np.float64]:
    reference = _as_binary_matrix(reference_matrix)
    if reference.shape[0] == 0:
        raise ValueError("reference_matrix must contain at least one individual.")
    saleem_ref = reference.T.astype(np.float64, copy=False)
    return (saleem_ref @ saleem_ref.T) / float(saleem_ref.shape[1])


def _sokal_michener(matrix: Any) -> Any:
    n_individuals = matrix.shape[1]
    if n_individuals == 0:
        raise ValueError("matrix must contain at least one individual column.")
    return matrix @ matrix.T / float(n_individuals)


def _resolve_torch_device(torch_module: Any, device: str | None) -> Any:
    if device is not None:
        return torch_module.device(device)
    if torch_module.cuda.is_available():
        return torch_module.device("cuda")
    if getattr(torch_module.backends, "mps", None) and torch_module.backends.mps.is_available():
        return torch_module.device("mps")
    return torch_module.device("cpu")


def _run_adam_epochs(
    tensor: Any,
    *,
    steps: int,
    learning_rate: float,
    loss_fn: Any,
    torch_module: Any,
) -> tuple[Any, float]:
    if steps <= 0:
        loss_value = float(loss_fn(tensor).detach().cpu().item())
        return tensor, loss_value
    variable = tensor.detach().requires_grad_(True)
    optimizer = torch_module.optim.Adam([variable], lr=learning_rate)
    loss_value = float("nan")
    for _ in range(steps):
        optimizer.zero_grad(set_to_none=True)
        loss = loss_fn(variable)
        loss.backward()
        optimizer.step()
        loss_value = float(loss.detach().cpu().item())
    return variable.detach(), loss_value


def _frequency_loss(reconstructed: Any, target_frequencies: Any) -> Any:
    current = reconstructed.sum(dim=1)
    return ((current - target_frequencies) ** 2).mean()


def _corr_loss_fn(target_corr: Any, torch_module: Any) -> Any:
    def loss(reconstructed: Any) -> Any:
        current_corr = _sokal_michener(reconstructed)
        return torch_module.norm(current_corr - target_corr, p="fro")

    return loss


def _count_preserving_binarize(
    continuous: ArrayLike,
    target_counts: NDArray[np.integer[Any]],
    *,
    rng: np.random.Generator,
) -> NDArray[np.int8]:
    matrix = np.asarray(continuous, dtype=np.float64)
    if matrix.ndim != 2:
        raise ValueError("continuous must be a two-dimensional matrix.")
    if matrix.shape[0] != target_counts.shape[0]:
        raise ValueError("target_counts length must match continuous row count.")

    binary = np.zeros(matrix.shape, dtype=np.int8)
    column_count = matrix.shape[1]
    tie_breakers = rng.random(matrix.shape)
    scores = matrix + tie_breakers * 1e-12
    for row_idx, count in enumerate(target_counts.astype(int)):
        carrier_count = int(count)
        if carrier_count <= 0:
            continue
        if carrier_count >= column_count:
            binary[row_idx, :] = 1
            continue
        top_columns = np.argpartition(-scores[row_idx], carrier_count - 1)[:carrier_count]
        binary[row_idx, top_columns] = 1
    return binary


def _binarize_reconstructed(
    continuous: ArrayLike,
    target_counts: NDArray[np.integer[Any]],
    *,
    rng: np.random.Generator,
) -> NDArray[np.int8]:
    return _count_preserving_binarize(continuous, target_counts, rng=rng)


def _query_initialization(
    target_matrix: ArrayLike,
    *,
    rng: np.random.Generator,
) -> NDArray[np.int8]:
    matrix = _as_binary_matrix(target_matrix, "target_matrix")
    reconstructed = np.zeros_like(matrix)
    column_count = matrix.shape[1]
    for row_idx, row in enumerate(matrix):
        carrier_count = int(np.sum(row))
        if carrier_count <= 0:
            continue
        columns = rng.choice(column_count, size=carrier_count, replace=False)
        reconstructed[row_idx, columns] = 1
    return reconstructed


def _compare_and_sort_columns(
    reference_matrix: ArrayLike,
    candidate_matrix: ArrayLike,
) -> NDArray[np.int8]:
    from scipy.spatial.distance import hamming

    reference = _as_binary_matrix(reference_matrix, "reference_matrix")
    candidate = _as_binary_matrix(candidate_matrix, "candidate_matrix")
    if reference.shape != candidate.shape:
        raise ValueError("reference_matrix and candidate_matrix must have the same shape.")
    aligned = np.zeros_like(reference)
    for column_idx in range(reference.shape[1]):
        distances = [
            (hamming(reference[:, column_idx].ravel(), candidate[:, other_idx].ravel()), other_idx)
            for other_idx in range(candidate.shape[1])
        ]
        _, best_idx = min(distances)
        aligned[:, column_idx] = candidate[:, best_idx]
    return aligned


def _official_repo_initialization(
    ground_truth_beacon: ArrayLike,
    *,
    rng: np.random.Generator,
) -> NDArray[np.int8]:
    true_saleem = _to_saleem_orientation(ground_truth_beacon)
    queried = _query_initialization(true_saleem, rng=rng)
    aligned = _compare_and_sort_columns(true_saleem, queried)
    return _from_saleem_orientation(aligned)


def _run_saleem_optimization(
    target_counts: NDArray[np.integer[Any]],
    reference_correlation: NDArray[np.float64],
    *,
    beacon_size: int,
    corr_epoch: int,
    freq_epoch: int,
    max_iterations: int,
    flip_threshold: int,
    min_outer_iterations: int,
    learning_rate: float,
    rng: np.random.Generator,
    torch_module: Any,
    device: Any,
    verbose: bool,
    initial_matrix: NDArray[np.int8],
) -> tuple[NDArray[np.int8], dict[str, Any]]:
    variant_count = int(target_counts.shape[0])
    if reference_correlation.shape != (variant_count, variant_count):
        raise ValueError("reference_correlation shape must match variant count.")

    target_frequencies = torch_module.tensor(
        target_counts.astype(np.float32), dtype=torch_module.float32, device=device
    )
    target_corr = torch_module.tensor(
        reference_correlation.astype(np.float32), dtype=torch_module.float32, device=device
    )
    corr_loss_fn = _corr_loss_fn(target_corr, torch_module)
    freq_loss_fn = lambda tensor: _frequency_loss(tensor, target_frequencies)

    initial = _as_binary_matrix(initial_matrix, "initial_matrix")
    if initial.shape != (variant_count, beacon_size):
        raise ValueError("initial_matrix shape must be (variant_count, beacon_size).")
    reconstructed = initial.copy()
    prev_reconstructed: NDArray[np.int8] | None = None
    outer_iterations = 0
    converged = False
    corr_loss_last = float("nan")
    freq_loss_last = float("nan")
    reconstructed_tensor = torch_module.tensor(
        reconstructed.astype(np.float32), dtype=torch_module.float32, device=device
    )

    for outer_iterations in range(1, max_iterations + 1):
        if verbose:
            print(
                f"saleem: outer_iter={outer_iterations}/{max_iterations} "
                f"(corr_epoch={corr_epoch}, freq_epoch={freq_epoch}, device={device})",
                flush=True,
            )

        reconstructed_tensor, corr_loss_last = _run_adam_epochs(
            reconstructed_tensor,
            steps=corr_epoch,
            learning_rate=learning_rate,
            loss_fn=corr_loss_fn,
            torch_module=torch_module,
        )
        reconstructed_tensor, freq_loss_last = _run_adam_epochs(
            reconstructed_tensor,
            steps=freq_epoch,
            learning_rate=learning_rate,
            loss_fn=freq_loss_fn,
            torch_module=torch_module,
        )

        reconstructed = _binarize_reconstructed(
            reconstructed_tensor.detach().cpu().numpy(),
            target_counts,
            rng=rng,
        )
        reconstructed_tensor = torch_module.tensor(
            reconstructed.astype(np.float32), dtype=torch_module.float32, device=device
        )

        if prev_reconstructed is not None and outer_iterations >= min_outer_iterations:
            flips = int(np.sum(reconstructed != prev_reconstructed))
            if flips < flip_threshold:
                converged = True
                break
        prev_reconstructed = reconstructed.copy()

    metadata = {
        "outer_iterations": int(outer_iterations),
        "converged": bool(converged),
        "corr_epoch": int(corr_epoch),
        "freq_epoch": int(freq_epoch),
        "max_iterations": int(max_iterations),
        "flip_threshold": int(flip_threshold),
        "min_outer_iterations": int(min_outer_iterations),
        "learning_rate": float(learning_rate),
        "final_corr_loss": corr_loss_last,
        "final_freq_loss": freq_loss_last,
        "device": str(device),
    }
    return reconstructed, metadata


def saleem_reconstruction_attack(
    problem: BeaconAttackProblem,
    reference_matrix: ArrayLike,
    *,
    prior_id: str = "saleem_ld_reference",
    initialization_seed: int | None = None,
    optimizer_seed: int | None = None,
    ground_truth_beacon: ArrayLike,
    corr_epoch: int = 1000,
    freq_epoch: int = 500,
    max_iterations: int = 3000,
    flip_threshold: int = 10,
    min_outer_iterations: int = 1,
    learning_rate: float = 0.001,
    device: str | None = None,
    verbose: bool = False,
) -> ReconstructionResult:
    """Run Saleem et al. (2025) alternating Adam reconstruction on carrier counts."""

    torch_module = _require_torch()
    attack_problem = _ensure_attack_problem(problem, "saleem_reconstruction_attack")
    if corr_epoch < 0 or freq_epoch < 0 or max_iterations < 1 or flip_threshold < 0:
        raise ValueError("Saleem optimization settings must be non-negative / positive as required.")
    if not np.isfinite(learning_rate) or learning_rate <= 0:
        raise ValueError("learning_rate must be a positive finite value.")

    reference = _as_binary_matrix(reference_matrix, "reference_matrix")
    variant_count = len(attack_problem.variant_ids)
    if reference.shape[1] != variant_count:
        raise ValueError("reference_matrix column count must match problem variants.")

    truth = _as_binary_matrix(ground_truth_beacon, "ground_truth_beacon")
    if truth.shape != (attack_problem.beacon_size, variant_count):
        raise ValueError("ground_truth_beacon shape must match beacon_size × variant_count.")

    target_counts = np.asarray(attack_problem.target_counts, dtype=np.int64)
    init_seed = _resolve_seed(initialization_seed, attack_problem.config.initialization_seed)
    opt_seed = _resolve_seed(optimizer_seed, attack_problem.config.optimizer_seed)
    rng = np.random.default_rng(init_seed)
    torch_module.manual_seed(opt_seed)
    resolved_device = _resolve_torch_device(torch_module, device)

    oracle_init = _official_repo_initialization(truth, rng=rng)
    initial_matrix = _to_saleem_orientation(oracle_init)
    reference_correlation = _reference_panel_to_sokal_michener(reference)
    saleem_reconstructed, optimization_metadata = _run_saleem_optimization(
        target_counts,
        reference_correlation,
        beacon_size=attack_problem.beacon_size,
        corr_epoch=corr_epoch,
        freq_epoch=freq_epoch,
        max_iterations=max_iterations,
        flip_threshold=flip_threshold,
        min_outer_iterations=min_outer_iterations,
        learning_rate=learning_rate,
        rng=rng,
        torch_module=torch_module,
        device=resolved_device,
        verbose=verbose,
        initial_matrix=initial_matrix,
    )
    reconstructed = _from_saleem_orientation(saleem_reconstructed)
    observed_counts = reconstructed.sum(axis=0, dtype=np.int64)
    max_count_error = int(np.max(np.abs(observed_counts - target_counts), initial=0))

    seeds = {"initialization_seed": int(init_seed), "optimizer_seed": int(opt_seed)}
    if attack_problem.config.group_selection_seed is not None:
        seeds["group_selection_seed"] = int(attack_problem.config.group_selection_seed)

    metadata: dict[str, Any] = {
        "method_name": SALEEM_RECONSTRUCTION_ATTACK,
        "manifest_id": attack_problem.manifest_id,
        "split_id": attack_problem.split_id,
        "prior_id": prior_id,
        "source_type": "saleem_official_algorithm",
        "seeds": dict(seeds),
        "max_count_error": max_count_error,
        "initialization_source": "official_repo_query_sort",
        "reference_panel_rows": int(reference.shape[0]),
        "config": {
            "corr_epoch": int(corr_epoch),
            "freq_epoch": int(freq_epoch),
            "max_iterations": int(max_iterations),
            "flip_threshold": int(flip_threshold),
            "min_outer_iterations": int(min_outer_iterations),
            "learning_rate": float(learning_rate),
            "device": str(resolved_device),
        },
        "optimization": optimization_metadata,
    }

    return ReconstructionResult(
        manifest_id=attack_problem.manifest_id,
        split_id=attack_problem.split_id,
        prior_id=prior_id,
        method_name=SALEEM_RECONSTRUCTION_ATTACK,
        seeds=seeds,
        reconstructed_binary=reconstructed,
        metrics={"max_count_error": max_count_error},
        metadata=metadata,
    )
