"""Reference-panel helpers for control conditions."""

from __future__ import annotations

import numpy as np


def shuffle_reference_columns(reference_matrix: np.ndarray, *, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    shuffled = reference_matrix.copy()
    for column in range(shuffled.shape[1]):
        shuffled[:, column] = rng.permutation(shuffled[:, column])
    return shuffled


def sample_reference_ids(reference_ids: list[str], *, panel_size: int, seed: int) -> list[str]:
    rng = np.random.default_rng(seed)
    if panel_size > len(reference_ids):
        raise ValueError("panel_size exceeds available reference rows.")
    indices = rng.choice(len(reference_ids), size=panel_size, replace=False)
    return [reference_ids[int(index)] for index in indices]
