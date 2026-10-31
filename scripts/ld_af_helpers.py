"""Shared Saleem baseline helpers for Genome Biology experiment reruns."""

from __future__ import annotations

import argparse
import csv
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Any

from beacon_reconstruction.core import BeaconAttackProblem, ReconstructionResult
from beacon_reconstruction.ld_af_baseline import SALEEM_RECONSTRUCTION_ATTACK, saleem_reconstruction_attack


def add_saleem_arguments(parser: ArgumentParser) -> None:
    group = parser.add_argument_group("Saleem baseline")
    group.add_argument("--saleem-corr-epoch", type=int, default=1000)
    group.add_argument("--saleem-freq-epoch", type=int, default=500)
    group.add_argument("--saleem-max-iterations", type=int, default=3000)
    group.add_argument("--saleem-flip-threshold", type=int, default=10)
    group.add_argument("--saleem-min-outer-iterations", type=int, default=1)
    group.add_argument("--saleem-learning-rate", type=float, default=0.001)
    group.add_argument(
        "--run-saleem",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run Saleem reconstruction attack (default: true).",
    )
    group.add_argument(
        "--saleem-only",
        action="store_true",
        help="Run only Saleem; merge other methods from --merge-metrics.",
    )
    group.add_argument(
        "--merge-metrics",
        type=Path,
        default=None,
        help="Existing metrics CSV whose non-Saleem rows are preserved when using --saleem-only.",
    )


def saleem_config_dict(args: Namespace) -> dict[str, Any]:
    return {
        "saleem_corr_epoch": args.saleem_corr_epoch,
        "saleem_freq_epoch": args.saleem_freq_epoch,
        "saleem_max_iterations": args.saleem_max_iterations,
        "saleem_flip_threshold": args.saleem_flip_threshold,
        "saleem_min_outer_iterations": args.saleem_min_outer_iterations,
        "saleem_learning_rate": args.saleem_learning_rate,
        "baseline_method": SALEEM_RECONSTRUCTION_ATTACK,
    }


def run_saleem_attack(
    problem: BeaconAttackProblem,
    public_ref_matrix: Any,
    beacon_matrix: Any,
    args: Namespace,
) -> ReconstructionResult:
    return saleem_reconstruction_attack(
        problem,
        public_ref_matrix,
        ground_truth_beacon=beacon_matrix,
        corr_epoch=args.saleem_corr_epoch,
        freq_epoch=args.saleem_freq_epoch,
        max_iterations=args.saleem_max_iterations,
        flip_threshold=args.saleem_flip_threshold,
        min_outer_iterations=args.saleem_min_outer_iterations,
        learning_rate=args.saleem_learning_rate,
    )


def read_metrics_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def merge_saleem_rows(
    *,
    new_rows: list[dict[str, Any]],
    merge_metrics: Path | None,
    saleem_only: bool,
    method_key: str = "method",
    saleem_method: str = SALEEM_RECONSTRUCTION_ATTACK,
    saleem_methods: tuple[str, ...] | None = None,
) -> list[dict[str, Any]]:
    if not saleem_only:
        return new_rows
    if merge_metrics is None:
        raise ValueError("--merge-metrics is required when --saleem-only is set.")
    if not merge_metrics.exists():
        raise FileNotFoundError(f"Missing merge metrics file: {merge_metrics}")
    methods_to_replace = saleem_methods or (saleem_method,)
    kept = [
        row for row in read_metrics_csv(merge_metrics) if row[method_key] not in methods_to_replace
    ]
    return kept + new_rows


def merge_saleem_overlap_rows(
    *,
    new_rows: list[dict[str, Any]],
    merge_metrics: Path | None,
    saleem_only: bool,
    overlap_filename: str,
    method_key: str = "method",
    saleem_method: str = SALEEM_RECONSTRUCTION_ATTACK,
    saleem_methods: tuple[str, ...] | None = None,
) -> list[dict[str, Any]]:
    if not saleem_only or merge_metrics is None:
        return new_rows
    overlap_path = merge_metrics.parent / overlap_filename
    if not overlap_path.exists():
        return new_rows
    methods_to_replace = saleem_methods or (saleem_method,)
    kept = [
        row for row in read_metrics_csv(overlap_path) if row[method_key] not in methods_to_replace
    ]
    return kept + new_rows
