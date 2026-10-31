"""Compare Stage 1 greedy templates with the exact minimum-L1 template.

R and N are small enough to enumerate every multiset of N reference rows.
The optimum is

    T* = argmin_T ||y - c(T)||_1

over that same with-replacement family. Reported quantities are the absolute
gap, the relative gap, and the optimal-hit rate.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from beacon_reconstruction.reconstruction import _select_reference_templates, minimum_template_count_l1
from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers, matrix_for_samples as _matrix_for_samples
from scripts.lib.sample_splits import load_sample_manifest, load_seed_manifest, split_samples


MANIFEST_ID = "1000g_eur_hg38_greedy_template_optimality_v1"
DEFAULT_MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
DEFAULT_OUT_DIR = ROOT / "results" / "block_inheritance_v1" / "greedy_template_optimality"
DEFAULT_SETTINGS = "5:5,6:4,8:5,8:8,10:6,10:8"
MAX_SEARCH_SPACE = 100_000
BOOTSTRAP_SEED = 0
BOOTSTRAP_ITERATIONS = 10_000


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}.")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _parse_settings(value: str) -> list[tuple[int, int]]:
    settings: list[tuple[int, int]] = []
    for item in value.split(","):
        text = item.strip()
        if not text:
            continue
        reference_size, beacon_size = text.split(":")
        settings.append((int(reference_size), int(beacon_size)))
    if not settings:
        raise ValueError("At least one R:N setting is required.")
    return settings


def _search_space_size(reference_size: int, beacon_size: int) -> int:
    return math.comb(reference_size + beacon_size - 1, beacon_size)


def _gaps(greedy_l1: int, optimal_l1: int) -> tuple[int, float, int]:
    absolute_gap = int(greedy_l1) - int(optimal_l1)
    relative_gap = float(absolute_gap) / float(max(1, int(optimal_l1)))
    optimal_hit = int(absolute_gap == 0)
    return absolute_gap, relative_gap, optimal_hit


def _mean(values: list[float]) -> float:
    return float(np.mean(np.asarray(values, dtype=float)))


def _bootstrap_mean_ci(values: list[float], *, seed: int) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    if array.size == 0:
        return float("nan"), float("nan")
    if array.size == 1:
        value = float(array[0])
        return value, value
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, array.size, size=(BOOTSTRAP_ITERATIONS, array.size))
    means = array[indices].mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def _aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    keys = sorted({(int(row["reference_size"]), int(row["beacon_size"])) for row in rows})
    for reference_size, beacon_size in keys:
        selected = [
            row
            for row in rows
            if int(row["reference_size"]) == reference_size and int(row["beacon_size"]) == beacon_size
        ]
        absolute = [float(row["absolute_gap"]) for row in selected]
        relative = [float(row["relative_gap"]) for row in selected]
        hits = [float(row["optimal_hit"]) for row in selected]
        hit_low, hit_high = _bootstrap_mean_ci(hits, seed=BOOTSTRAP_SEED + reference_size * 100 + beacon_size)
        out.append(
            {
                "reference_size": reference_size,
                "beacon_size": beacon_size,
                "search_space": selected[0]["search_space"],
                "n_instances": len(selected),
                "mean_greedy_l1": _mean([float(row["greedy_l1"]) for row in selected]),
                "mean_optimal_l1": _mean([float(row["optimal_l1"]) for row in selected]),
                "mean_absolute_gap": _mean(absolute),
                "median_absolute_gap": float(np.median(absolute)),
                "max_absolute_gap": float(np.max(absolute)),
                "mean_relative_gap": _mean(relative),
                "median_relative_gap": float(np.median(relative)),
                "optimal_hit_rate": _mean(hits),
                "optimal_hit_ci95_low": hit_low,
                "optimal_hit_ci95_high": hit_high,
            }
        )
    return out


def _write_report(path: Path, aggregate_rows: list[dict[str, Any]], config: dict[str, Any]) -> None:
    lines = [
        "# Exact Template Optimality At Small R And N",
        "",
        "## Configuration",
        "",
        f"- windows: {config['max_windows']}",
        f"- replicates: {config['max_replicates']}",
        f"- draws per window and replicate: {config['draws']}",
        f"- settings R:N: {config['settings']}",
        "- reference panel: shared matched panel, disjoint from the beacon",
        "- beacon rows: random subset of the N=25 cohort",
        "- feasible templates: multisets of N rows from R reference rows, with replacement",
        "",
        "`T*` minimizes `||y - c(T)||_1` over that family.",
        "Absolute gap is `||y - c(T_greedy)||_1 - ||y - c(T*)||_1`.",
        "Relative gap divides the absolute gap by `max(1, ||y - c(T*)||_1)`.",
        "Optimal-hit rate is the fraction of instances with absolute gap 0.",
        f"The hit-rate interval is a percentile bootstrap of that fraction ({BOOTSTRAP_ITERATIONS} resamples, seed {BOOTSTRAP_SEED}).",
        "",
        "## Aggregate",
        "",
        "| R | N | space | n | mean greedy L1 | mean optimal L1 | mean abs gap | median abs gap | max abs gap | mean rel gap | hit rate | hit 95% CI |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['reference_size']} | {row['beacon_size']} | {row['search_space']} | {row['n_instances']} | "
            f"{row['mean_greedy_l1']:.2f} | {row['mean_optimal_l1']:.2f} | {row['mean_absolute_gap']:.2f} | "
            f"{row['median_absolute_gap']:.2f} | {row['max_absolute_gap']:.0f} | {row['mean_relative_gap']:.4f} | "
            f"{row['optimal_hit_rate']:.4f} | [{row['optimal_hit_ci95_low']:.4f}, {row['optimal_hit_ci95_high']:.4f}] |"
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "This comparison is restricted to enumerated small panels. It evaluates the final count objective",
            "of Stage 1 template selection and does not include the Stage 2 projection.",
            "",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--draws", type=int, default=5)
    parser.add_argument("--settings", default=DEFAULT_SETTINGS, help="Comma-separated R:N pairs.")
    parser.add_argument("--multiwindow-dir", type=Path, default=DEFAULT_MULTIWINDOW_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    if args.draws < 1:
        parser.error("--draws must be positive")

    settings = _parse_settings(args.settings)
    for reference_size, beacon_size in settings:
        if reference_size < 1 or beacon_size < 1:
            parser.error(f"R and N must be positive, got {reference_size}:{beacon_size}")
        if beacon_size > 25:
            parser.error("N must be at most 25 so the beacon subset fits in the N=25 cohort.")
        space = _search_space_size(reference_size, beacon_size)
        if space > MAX_SEARCH_SPACE:
            parser.error(f"R={reference_size}, N={beacon_size} has {space} multisets, above {MAX_SEARCH_SPACE}.")

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    windows = _read_csv(args.multiwindow_dir / "multiwindow_manifest.csv")[: max(0, args.max_windows)]
    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]
    config = {
        "manifest_id": MANIFEST_ID,
        "max_windows": args.max_windows,
        "max_replicates": args.max_replicates,
        "draws": args.draws,
        "settings": args.settings,
        "max_search_space": MAX_SEARCH_SPACE,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
        "reference_panel_mode": "shared",
    }

    rows: list[dict[str, Any]] = []
    for window in windows:
        window_id = window["window_id"]
        sample_ids, _variant_ids, _positions, _target_alleles, carrier_matrix = _load_raw_as_target_carriers(
            Path(window["raw_path"])
        )
        for replicate_id in replicate_ids:
            seeds = seeds_by_replicate[replicate_id]
            split = split_samples(sample_rows, replicate_id, 25, shared_reference_panel=True)
            beacon_ids = split["beacon"]
            reference_ids = split["public_ref"]
            for reference_size, beacon_size in settings:
                space = _search_space_size(reference_size, beacon_size)
                for draw in range(args.draws):
                    subset_rng = np.random.default_rng(seeds["group_selection_seed"] + 10007 * draw + beacon_size + 17 * reference_size)
                    beacon_positions = subset_rng.choice(len(beacon_ids), size=beacon_size, replace=False)
                    reference_positions = subset_rng.choice(len(reference_ids), size=reference_size, replace=False)
                    beacon_subset = [beacon_ids[int(index)] for index in beacon_positions]
                    reference_subset = [reference_ids[int(index)] for index in reference_positions]
                    beacon_matrix = _matrix_for_samples(sample_ids, carrier_matrix, beacon_subset)
                    reference_matrix = _matrix_for_samples(sample_ids, carrier_matrix, reference_subset)
                    target_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                    greedy = _select_reference_templates(
                        reference_matrix,
                        target_counts,
                        beacon_size=beacon_size,
                        rng=np.random.default_rng(seeds["initialization_seed"] + draw),
                    )
                    greedy_l1 = int(np.abs(greedy.sum(axis=0, dtype=np.int64) - target_counts).sum())
                    optimal_l1 = minimum_template_count_l1(
                        reference_matrix,
                        target_counts,
                        beacon_size=beacon_size,
                    )
                    absolute_gap, relative_gap, optimal_hit = _gaps(greedy_l1, optimal_l1)
                    if absolute_gap < 0:
                        raise RuntimeError(f"Greedy beat the enumerated optimum for {window_id} {replicate_id}.")
                    rows.append(
                        {
                            "window_id": window_id,
                            "replicate_id": replicate_id,
                            "draw": draw,
                            "reference_size": reference_size,
                            "beacon_size": beacon_size,
                            "search_space": space,
                            "variant_count": int(target_counts.shape[0]),
                            "greedy_l1": greedy_l1,
                            "optimal_l1": optimal_l1,
                            "absolute_gap": absolute_gap,
                            "relative_gap": relative_gap,
                            "optimal_hit": optimal_hit,
                            "source_overlap_count": len(set(beacon_subset).intersection(reference_subset)),
                        }
                    )

    aggregate_rows = _aggregate(rows)
    _write_csv(out_dir / "greedy_template_optimality_metrics.csv", rows)
    _write_csv(out_dir / "greedy_template_optimality_aggregate.csv", aggregate_rows)
    _write_json(out_dir / "greedy_template_optimality_manifest.json", config)
    _write_report(out_dir / "greedy_template_optimality_report.md", aggregate_rows, config)


if __name__ == "__main__":
    main()
