"""Load multi-window PLINK carrier matrices for Beacon experiments."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def _raw_column_key(column: str) -> tuple[str, str]:
    variant_id, allele = column.rsplit("_", 1)
    return variant_id, allele


def _other_allele(variant_id: str, allele: str) -> str:
    parts = variant_id.split(":")
    if len(parts) < 4:
        raise ValueError(f"Cannot infer alleles from variant_id={variant_id}.")
    ref, alt = parts[-2], parts[-1]
    if allele == ref:
        return alt
    if allele == alt:
        return ref
    raise ValueError(f"Raw allele {allele} is not present in variant_id={variant_id}.")


def load_raw_carriers(
    raw_path: Path,
    *,
    target_alleles: dict[str, str] | None = None,
) -> tuple[list[str], tuple[str, ...], tuple[int, ...], dict[str, str], np.ndarray]:
    with raw_path.open() as handle:
        header = handle.readline().strip().split()
        raw_columns = header[6:]
        parsed = [_raw_column_key(column) for column in raw_columns]
        sample_ids: list[str] = []
        dosage_rows: list[np.ndarray] = []
        for line in handle:
            parts = line.strip().split()
            sample_ids.append(parts[1])
            dosage_rows.append(np.asarray(parts[6:], dtype=float))

    dosage_matrix = np.vstack(dosage_rows)
    all_variant_ids = tuple(variant_id for variant_id, _allele in parsed)
    if target_alleles is None:
        target_alleles = {}
        selected_columns = list(range(len(raw_columns)))
        for source_index, (variant_id, raw_allele) in enumerate(parsed):
            raw_carrier_count = int(np.sum(dosage_matrix[:, source_index] > 0))
            other_carrier_count = int(np.sum(dosage_matrix[:, source_index] < 2))
            target_alleles[variant_id] = (
                raw_allele if raw_carrier_count <= other_carrier_count else _other_allele(variant_id, raw_allele)
            )
    else:
        selected_columns = [all_variant_ids.index(variant_id) for variant_id in target_alleles]

    variant_ids = tuple(target_alleles.keys())
    positions = tuple(int(variant_id.split(":")[1]) for variant_id in variant_ids)
    carriers = np.zeros((len(sample_ids), len(variant_ids)), dtype=np.int8)
    for out_index, source_index in enumerate(selected_columns):
        variant_id, raw_allele = parsed[source_index]
        target_allele = target_alleles[variant_id]
        if raw_allele == target_allele:
            carriers[:, out_index] = dosage_matrix[:, source_index] > 0
        else:
            carriers[:, out_index] = dosage_matrix[:, source_index] < 2
    return sample_ids, variant_ids, positions, dict(target_alleles), carriers.astype(np.int8)


def matrix_for_samples(sample_ids: list[str], carrier_matrix: np.ndarray, selected_ids: list[str]) -> np.ndarray:
    index = {sample_id: i for i, sample_id in enumerate(sample_ids)}
    return carrier_matrix[[index[sample_id] for sample_id in selected_ids], :]
