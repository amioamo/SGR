"""Sample manifests and beacon/reference splits for 1000 Genomes experiments."""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST_DIR = ROOT / "data" / "manifests"
MANIFEST_ID = "1000g_eur_hg38_chr22_dense500_v1"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def load_sample_manifest() -> list[dict[str, str]]:
    rows = _read_csv(MANIFEST_DIR / "1000g_eur_chr22_prototype_samples.csv")
    if {row["manifest_id"] for row in rows} != {MANIFEST_ID}:
        raise ValueError("Sample manifest does not match expected manifest ID.")
    return rows


def load_seed_manifest() -> dict[str, dict[str, int]]:
    rows = _read_csv(MANIFEST_DIR / "1000g_eur_chr22_prototype_seeds.csv")
    seeds: dict[str, dict[str, int]] = {}
    for row in rows:
        seeds[row["replicate_id"]] = {
            "cohort_seed": int(row["cohort_seed"]),
            "variant_window_seed": int(row["variant_window_seed"]),
            "initialization_seed": int(row["initialization_seed"]),
            "optimizer_seed": int(row["optimizer_seed"]),
            "group_selection_seed": int(row["group_selection_seed"]),
        }
    return seeds


def split_samples(
    sample_rows: list[dict[str, str]],
    replicate_id: str,
    beacon_size: int,
    *,
    shared_reference_panel: bool = False,
) -> dict[str, list[str]]:
    rows = [row for row in sample_rows if row["replicate_id"] == replicate_id]
    beacon_n25 = [
        row["sample_id"]
        for row in rows
        if row["role"] == "beacon_cohort" and row["nested_beacon_group"] == "beacon_N25"
    ]
    beacon_n50_only = [
        row["sample_id"]
        for row in rows
        if row["role"] == "beacon_cohort" and row["nested_beacon_group"] == "beacon_N50_only"
    ]
    beacon_n100_only = [
        row["sample_id"]
        for row in rows
        if row["role"] == "beacon_cohort" and row["nested_beacon_group"] == "beacon_N100_only"
    ]
    beacon_n250_only = [
        row["sample_id"]
        for row in rows
        if row["role"] == "beacon_cohort" and row["nested_beacon_group"] == "beacon_N250_only"
    ]
    if beacon_size == 25:
        beacon = beacon_n25
    elif beacon_size == 50:
        beacon = beacon_n25 + beacon_n50_only
    elif beacon_size == 100:
        beacon = beacon_n25 + beacon_n50_only + beacon_n100_only
    elif beacon_size == 250:
        beacon = beacon_n25 + beacon_n50_only + beacon_n100_only + beacon_n250_only
    else:
        raise ValueError("Only beacon sizes 25, 50, 100, and 250 are supported.")
    public_ref = [row["sample_id"] for row in rows if row["role"] == "public_af_ld_reference"]
    reference_panel = [row["sample_id"] for row in rows if row["role"] == "arg_reference_or_proxy_reference"]

    if len(beacon) != beacon_size or len(public_ref) != 100 or len(reference_panel) != 100:
        raise ValueError(f"Unexpected split sizes for replicate {replicate_id}, N={beacon_size}.")
    if set(beacon) & set(public_ref) or set(beacon) & set(reference_panel):
        raise ValueError(f"Source overlap detected for replicate {replicate_id}, N={beacon_size}.")
    if shared_reference_panel:
        reference_panel = list(public_ref)
    elif set(public_ref) & set(reference_panel):
        raise ValueError(f"Source overlap detected for replicate {replicate_id}, N={beacon_size}.")
    return {"beacon": beacon, "public_ref": public_ref, "arg_ref": reference_panel}
