"""Supplementary anchor: random |M'|=1000 SNP sets in the Saleem benchmark regime.

Pools variants from the six dense500 windows (3000 SNPs total), draws independent
random 1000-SNP subsets, and compares count-only / local LD refinement / Saleem LD optimization /
matched template-count synthesis under the shared R=100 reference panel.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, pstdev
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from beacon_reconstruction.baselines import independent_count, pairwise_moments, pairwise_refinement
from beacon_reconstruction.reconstruction import reference_haplotype_reconstruction
from beacon_reconstruction.core import ExperimentConfig, BeaconAttackProblem, BeaconEvaluationTruth, ReconstructionResult
from beacon_reconstruction.evaluation import EvaluationResult, evaluate_reconstruction
from beacon_reconstruction.ld_af_baseline import SALEEM_RECONSTRUCTION_ATTACK
from scripts.experiment_runtime_helpers import attach_reference_panel_mode, run_timed, runtime_metadata_fields
from scripts.lib.multiwindow_data import load_raw_carriers, matrix_for_samples
from scripts.lib.sample_splits import load_sample_manifest, load_seed_manifest, split_samples
from scripts.lib.reference_helpers import sample_reference_ids
from scripts.ld_af_helpers import add_saleem_arguments, run_saleem_attack, saleem_config_dict


MANIFEST_ID = "1000g_eur_hg38_saleem_random_m1000_anchor_v1"
DEFAULT_MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
DEFAULT_OUT_DIR = ROOT / "results" / "experiments" / "saleem_random_m1000_anchor"
POOL_WINDOW_ID = "pooled_dense500_x6"
DEFAULT_METHODS = (
    "independent_count",
    "pairwise_refinement",
    SALEEM_RECONSTRUCTION_ATTACK,
    "matched_template_n100",
)


@dataclass(frozen=True)
class PooledVariant:
    window_id: str
    variant_id: str
    position: int
    column_index: int

    @property
    def global_key(self) -> str:
        return f"{self.window_id}|{self.variant_id}"


@dataclass
class WindowData:
    window_id: str
    sample_ids: list[str]
    variant_ids: tuple[str, ...]
    positions: tuple[int, ...]
    carrier_matrix: np.ndarray


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _parse_methods(value: str) -> list[str]:
    methods = [item.strip() for item in value.split(",") if item.strip()]
    unknown = set(methods) - set(DEFAULT_METHODS)
    if unknown:
        raise ValueError(f"Unknown methods: {sorted(unknown)}; expected subset of {DEFAULT_METHODS}.")
    return methods


def _parse_ints(value: str) -> list[int]:
    return sorted({int(item.strip()) for item in value.split(",") if item.strip()})


def _serialize_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, dict):
        return {str(key): _serialize_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize_value(item) for item in value]
    return value


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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_serialize_value(payload), indent=2, sort_keys=True) + "\n")


def _load_windows(multiwindow_dir: Path, *, max_windows: int) -> list[WindowData]:
    windows = _read_csv(multiwindow_dir / "multiwindow_manifest.csv")[:max(0, max_windows)]
    loaded: list[WindowData] = []
    for window in windows:
        sample_ids, variant_ids, positions, _, carrier_matrix = load_raw_carriers(
            Path(window["raw_path"])
        )
        loaded.append(
            WindowData(
                window_id=window["window_id"],
                sample_ids=sample_ids,
                variant_ids=variant_ids,
                positions=positions,
                carrier_matrix=carrier_matrix.astype(np.int8, copy=False),
            )
        )
    reference_ids = loaded[0].sample_ids
    for item in loaded[1:]:
        if item.sample_ids != reference_ids:
            raise ValueError("All windows must share identical sample_id order for pooled SNP sets.")
    return loaded


def _build_pool(windows: list[WindowData]) -> list[PooledVariant]:
    pool: list[PooledVariant] = []
    for window in windows:
        for column_index, (variant_id, position) in enumerate(zip(window.variant_ids, window.positions)):
            pool.append(
                PooledVariant(
                    window_id=window.window_id,
                    variant_id=variant_id,
                    position=int(position),
                    column_index=column_index,
                )
            )
    return pool


def _window_lookup(windows: list[WindowData]) -> dict[str, WindowData]:
    return {window.window_id: window for window in windows}


def _draw_snp_set(
    pool: list[PooledVariant],
    *,
    snp_count: int,
    snp_set_replicate: int,
    base_seed: int,
) -> list[PooledVariant]:
    if snp_count > len(pool):
        raise ValueError(f"Requested {snp_count} SNPs but pooled pool has only {len(pool)}.")
    rng = np.random.default_rng(base_seed + snp_set_replicate * 10_007)
    indices = np.sort(rng.choice(len(pool), size=snp_count, replace=False).astype(int))
    return [pool[int(index)] for index in indices]


def _subset_matrix(sample_ids: list[str], windows_by_id: dict[str, WindowData], selected: list[PooledVariant], sample_subset: list[str]) -> tuple[tuple[str, ...], tuple[int, ...], np.ndarray]:
    parts: list[np.ndarray] = []
    variant_ids: list[str] = []
    positions: list[int] = []
    for item in selected:
        window = windows_by_id[item.window_id]
        column = window.carrier_matrix[:, item.column_index]
        index = {sample_id: idx for idx, sample_id in enumerate(window.sample_ids)}
        part = np.asarray([column[index[sample_id]] for sample_id in sample_subset], dtype=np.int8)
        parts.append(part[:, None])
        variant_ids.append(item.global_key)
        positions.append(item.position)
    matrix = np.hstack(parts)
    return tuple(variant_ids), tuple(positions), matrix


def _with_config(result: ReconstructionResult, *, config: dict[str, Any], optimizer_seed: int) -> ReconstructionResult:
    seeds = dict(result.seeds)
    seeds.setdefault("optimizer_seed", int(optimizer_seed))
    metadata = dict(result.metadata)
    metadata["config"] = config
    return ReconstructionResult(
        manifest_id=result.manifest_id,
        split_id=result.split_id,
        prior_id=result.prior_id,
        method_name=result.method_name,
        seeds=seeds,
        reconstructed_binary=result.reconstructed_binary,
        metrics=result.metrics,
        metadata=metadata,
    )


def _record(
    *,
    snp_count: int,
    snp_set_replicate: int,
    beacon_size: int,
    replicate_id: str,
    method: str,
    reference_population: str,
    template_design: str,
    panel_size: int | None,
    selected_keys: tuple[str, ...],
    result: ReconstructionResult,
    evaluation: EvaluationResult,
    source_overlap_count: int,
    runtime_fields: dict[str, Any] | None = None,
) -> dict[str, Any]:
    metrics = evaluation.metrics
    row = {
        "snp_count": snp_count,
        "snp_set_replicate": snp_set_replicate,
        "beacon_size": beacon_size,
        "replicate_id": replicate_id,
        "window_id": POOL_WINDOW_ID,
        "method": method,
        "method_name": result.method_name,
        "reference_population": reference_population,
        "template_design": template_design,
        "panel_size": panel_size,
        "variant_count": len(selected_keys),
        "selected_variant_ids": ";".join(selected_keys),
        "source_overlap_count": source_overlap_count,
        "precision": metrics.precision,
        "recall": metrics.recall,
        "f1": metrics.f1,
        "pairwise_error": metrics.pairwise_error,
        "alignment_total_hamming_cost": metrics.alignment_total_hamming_cost,
        "max_count_error": metrics.max_count_error,
        "invalid_run": evaluation.invalid_run.is_invalid,
        "invalid_reasons": ";".join(evaluation.invalid_run.reasons),
    }
    if runtime_fields:
        row.update(runtime_fields)
    return row


def _aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    keys = sorted({(int(row["beacon_size"]), row["method"]) for row in rows})
    for beacon_size, method in keys:
        selected = [row for row in rows if int(row["beacon_size"]) == beacon_size and row["method"] == method]
        f1s = [float(row["f1"]) for row in selected]
        entry: dict[str, Any] = {
            "beacon_size": beacon_size,
            "method": method,
            "reference_population": selected[0]["reference_population"],
            "template_design": selected[0]["template_design"],
            "panel_size": selected[0]["panel_size"],
            "snp_count": selected[0]["snp_count"],
            "n_runs": len(selected),
            "mean_f1": mean(f1s),
            "std_f1": pstdev(f1s) if len(f1s) > 1 else 0.0,
            "median_f1": float(np.median(f1s)),
            "total_source_overlap_count": int(sum(int(row["source_overlap_count"]) for row in selected)),
            "invalid_run_count": int(sum(_as_bool(row["invalid_run"]) for row in selected)),
        }
        if "wall_time_sec" in selected[0]:
            wall_times = [float(row["wall_time_sec"]) for row in selected]
            entry["mean_wall_time_sec"] = mean(wall_times)
            entry["std_wall_time_sec"] = pstdev(wall_times) if len(wall_times) > 1 else 0.0
        out.append(entry)
    return out


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, np.integer)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes"}
    return bool(value)


def _normalize_metric_row(row: dict[str, Any]) -> dict[str, Any]:
    normalized = dict(row)
    if "invalid_run" in normalized:
        normalized["invalid_run"] = _as_bool(normalized["invalid_run"])
    for key in ("beacon_size", "snp_count", "snp_set_replicate", "source_overlap_count", "panel_size"):
        if key in normalized and normalized[key] not in ("", None):
            normalized[key] = int(float(normalized[key]))
    for key in ("f1", "precision", "recall", "pairwise_error", "wall_time_sec"):
        if key in normalized and normalized[key] not in ("", None):
            normalized[key] = float(normalized[key])
    return normalized


def _read_existing_metrics(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return [_normalize_metric_row(row) for row in _read_csv(path)]


def _read_existing_manifest(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _merge_beacon_sizes(existing: list[int], new: list[int]) -> list[int]:
    return sorted(set(existing) | set(new))


def _write_report(path: Path, *, aggregate_rows: list[dict[str, Any]], config: dict[str, Any]) -> None:
    lines = [
        "# Saleem Random |M'|=1000 Supplementary Anchor",
        "",
        "## Purpose",
        "",
        "Supplementary anchor to the Saleem benchmark regime: genome-wide-style random SNP",
        "sampling (|M'|=1000) pooled from six dense500 windows, compared under the same",
        "shared R=100 disjoint reference panel as the main manuscript.",
        "",
        "## Configuration",
        "",
        f"- pooled windows: {config['max_windows']}",
        f"- |M'|: {config['snp_count']}",
        f"- SNP-set replicates: {config['snp_set_replicates']}",
        f"- cohort replicates: {config['max_replicates']}",
        f"- beacon sizes: {','.join(str(n) for n in config['beacon_sizes'])}",
        f"- methods: {', '.join(config['methods'])}",
        f"- reference panel: shared R=100 (disjoint from beacon only)",
        f"- subset seed base: {config['subset_seed']}",
        "",
        "## SNP set construction",
        "",
        "1. Load six pre-selected dense500 windows (~500 SNPs each; 3000 total).",
        "2. Pool all variants with keys `window_id|variant_id` (uniform pool; no LD/AF pre-filter).",
        "3. For each SNP-set replicate `k`, draw 1000 indices without replacement using",
        "   `numpy.random.Generator(seed=subset_seed + k * 10007).choice(len(pool), 1000, replace=False)`.",
        "4. Sort selected indices for deterministic column order; genotype columns are concatenated across windows.",
        "",
        "This matches Saleem's random |M'|=1000 design in cardinality and sampling mechanism, but the",
        "candidate pool is regional (six dense loci) rather than genome-wide HapMap SNPs.",
        "",
        "## Aggregate F1 (mean ± SD)",
        "",
        "| N | method | mean_f1 ± std_f1 | mean_wall_time_sec | n_runs | invalid |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        wall = row.get("mean_wall_time_sec")
        wall_text = f"{wall:.2f}" if wall is not None else "—"
        lines.append(
            f"| {row['beacon_size']} | {row['method']} | "
            f"{row['mean_f1']:.4f} ± {row['std_f1']:.4f} | {wall_text} | {row['n_runs']} | {row['invalid_run_count']} |"
        )
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            "- SNPs are drawn uniformly without replacement from the local 3000-SNP dense-window pool;",
            "  this proxies Saleem's random 1000-SNP design but remains **regional**, not genome-wide.",
            "- Intended for **Supplementary Fig/Table** + Methods sensitivity paragraph only.",
            "- Do not treat as the primary manuscript evaluation axis (main text keeps dense500 windows).",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--multiwindow-dir", type=Path, default=DEFAULT_MULTIWINDOW_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--max-replicates", type=int, default=3)
    parser.add_argument("--snp-set-replicates", type=int, default=3)
    parser.add_argument("--beacon-sizes", default="25,50,100,250")
    parser.add_argument(
        "--methods",
        default=",".join(DEFAULT_METHODS),
        help=f"Comma-separated subset of: {', '.join(DEFAULT_METHODS)}",
    )
    parser.add_argument("--snp-count", type=int, default=1000)
    parser.add_argument("--subset-seed", type=int, default=20260608)
    parser.add_argument("--template-panel-size", type=int, default=100)
    parser.add_argument("--pairwise-epochs", type=int, default=20)
    parser.add_argument("--pairwise-proposals-per-epoch", type=int, default=5000)
    parser.add_argument("--pairwise-early-stopping", type=int, default=5)
    parser.add_argument(
        "--append",
        action="store_true",
        help="Merge new beacon sizes or methods into existing metrics/manifest in --out-dir.",
    )
    add_saleem_arguments(parser)
    args = parser.parse_args()

    run_beacon_sizes = _parse_ints(args.beacon_sizes)
    run_methods = _parse_methods(args.methods)
    metrics_path = args.out_dir / "saleem_random_m1000_anchor_metrics.csv"
    manifest_path = args.out_dir / "saleem_random_m1000_anchor_manifest.json"
    existing_rows = _read_existing_metrics(metrics_path) if args.append else []
    existing_manifest = _read_existing_manifest(manifest_path) if args.append else None
    if args.append:
        prior_sizes = (
            existing_manifest["config"].get("beacon_sizes", [])
            if existing_manifest and "config" in existing_manifest
            else sorted({int(row["beacon_size"]) for row in existing_rows})
        )
        beacon_sizes = _merge_beacon_sizes(prior_sizes, run_beacon_sizes)
        prior_methods = (
            existing_manifest["config"].get("methods", [])
            if existing_manifest and "config" in existing_manifest
            else sorted({row["method"] for row in existing_rows})
        )
        methods = sorted(
            set(prior_methods) | set(run_methods),
            key=lambda item: DEFAULT_METHODS.index(item) if item in DEFAULT_METHODS else 99,
        )
        existing_rows = [
            row
            for row in existing_rows
            if row["method"] not in run_methods or int(row["beacon_size"]) not in run_beacon_sizes
        ]
    else:
        beacon_sizes = run_beacon_sizes
        methods = run_methods
    snp_set_replicates = list(range(max(0, args.snp_set_replicates)))
    windows = _load_windows(args.multiwindow_dir, max_windows=args.max_windows)
    pool = _build_pool(windows)
    windows_by_id = _window_lookup(windows)
    sample_rows = load_sample_manifest()
    seeds_by_replicate = load_seed_manifest()
    replicate_ids = sorted(seeds_by_replicate.keys())[: max(0, args.max_replicates)]

    config = {
        "max_windows": args.max_windows,
        "max_replicates": args.max_replicates,
        "snp_set_replicates": args.snp_set_replicates,
        "beacon_sizes": beacon_sizes,
        "snp_count": args.snp_count,
        "subset_seed": args.subset_seed,
        "template_panel_size": args.template_panel_size,
        "reference_panel_mode": "shared",
        "methods": methods,
        "pairwise_epochs": args.pairwise_epochs,
        "pairwise_proposals_per_epoch": args.pairwise_proposals_per_epoch,
        "pairwise_early_stopping": args.pairwise_early_stopping,
        **saleem_config_dict(args),
        "multiwindow_dir": str(args.multiwindow_dir),
        "out_dir": str(args.out_dir),
    }

    snp_set_manifest: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = list(existing_rows)

    if args.append and existing_manifest and "snp_sets" in existing_manifest:
        snp_set_manifest = list(existing_manifest["snp_sets"])
        if len(snp_set_manifest) != len(snp_set_replicates):
            raise ValueError("Existing manifest snp_sets count does not match --snp-set-replicates.")

    for snp_set_replicate in snp_set_replicates:
        if snp_set_manifest and snp_set_replicate < len(snp_set_manifest):
            selected_keys = tuple(snp_set_manifest[snp_set_replicate]["selected_variant_ids"])
            key_set = set(selected_keys)
            selected = [pool_item for pool_item in pool if pool_item.global_key in key_set]
            if len(selected) != args.snp_count:
                raise ValueError(
                    f"SNP-set replicate {snp_set_replicate}: expected {args.snp_count} variants, got {len(selected)}."
                )
        else:
            selected = _draw_snp_set(
                pool,
                snp_count=args.snp_count,
                snp_set_replicate=snp_set_replicate,
                base_seed=args.subset_seed,
            )
            selected_keys = tuple(item.global_key for item in selected)
            snp_set_manifest.append(
                {
                    "snp_set_replicate": snp_set_replicate,
                    "snp_count": args.snp_count,
                    "variant_count": len(selected),
                    "selected_variant_ids": list(selected_keys),
                }
            )

        for replicate_id in replicate_ids:
            seeds = seeds_by_replicate[replicate_id]
            for beacon_size in run_beacon_sizes:
                split = split_samples(
                    sample_rows,
                    replicate_id,
                    beacon_size,
                    shared_reference_panel=True,
                )
                template_ref_ids = split["public_ref"]
                panel_ids = sample_reference_ids(
                    template_ref_ids,
                    panel_size=args.template_panel_size,
                    seed=seeds["group_selection_seed"] + snp_set_replicate,
                )
                beacon_matrix = _subset_matrix(windows[0].sample_ids, windows_by_id, selected, split["beacon"])[2]
                public_ref_matrix = _subset_matrix(windows[0].sample_ids, windows_by_id, selected, split["public_ref"])[2]
                template_matrix = _subset_matrix(windows[0].sample_ids, windows_by_id, selected, panel_ids)[2]
                variant_ids, positions, _ = _subset_matrix(
                    windows[0].sample_ids, windows_by_id, selected, split["beacon"]
                )
                target_counts = beacon_matrix.sum(axis=0, dtype=np.int64)
                reference_pairwise = pairwise_moments(public_ref_matrix)
                split_id = (
                    f"m{args.snp_count}_snprep{snp_set_replicate}_rep{replicate_id}_N{beacon_size}"
                )
                runtime_config = {
                    "window_id": POOL_WINDOW_ID,
                    "replicate_id": replicate_id,
                    "snp_set_replicate": snp_set_replicate,
                    "beacon_size": beacon_size,
                    "snp_count": args.snp_count,
                    "config": config,
                }
                problem = BeaconAttackProblem(
                    manifest_id=MANIFEST_ID,
                    split_id=split_id,
                    beacon_size=beacon_size,
                    variant_ids=variant_ids,
                    positions=positions,
                    query_responses=target_counts,
                    target_counts=target_counts,
                    public_ld=reference_pairwise,
                    config=ExperimentConfig(
                        initialization_seed=seeds["initialization_seed"] + snp_set_replicate,
                        optimizer_seed=seeds["optimizer_seed"],
                        group_selection_seed=seeds["group_selection_seed"],
                        metadata={
                            "window_id": POOL_WINDOW_ID,
                            "experiment": "saleem_random_m1000_anchor",
                        },
                    ),
                )
                truth = BeaconEvaluationTruth(MANIFEST_ID, split_id, beacon_matrix, tuple(split["beacon"]))

                run_specs: list[dict[str, Any]] = [
                    {
                        "method": "independent_count",
                        "reference_population": "none",
                        "template_design": "count_only",
                        "panel_size": None,
                        "reference_ids": [],
                        "run": lambda problem=problem, reference_pairwise=reference_pairwise, runtime_config=runtime_config, seeds=seeds: _with_config(
                            independent_count(problem, reference_pairwise=reference_pairwise),
                            config=runtime_config,
                            optimizer_seed=seeds["optimizer_seed"],
                        ),
                    },
                    {
                        "method": "pairwise_refinement",
                        "reference_population": "eur",
                        "template_design": "pairwise_ld",
                        "panel_size": len(split["public_ref"]),
                        "reference_ids": split["public_ref"],
                        "run": lambda problem=problem, reference_pairwise=reference_pairwise, runtime_config=runtime_config, seeds=seeds: _with_config(
                            pairwise_refinement(
                                problem,
                                reference_pairwise,
                                epochs=args.pairwise_epochs,
                                proposals_per_epoch=args.pairwise_proposals_per_epoch,
                                early_stopping=args.pairwise_early_stopping,
                            ),
                            config=runtime_config,
                            optimizer_seed=seeds["optimizer_seed"],
                        ),
                    },
                    {
                        "method": SALEEM_RECONSTRUCTION_ATTACK,
                        "reference_population": "eur",
                        "template_design": "saleem_ld",
                        "panel_size": len(split["public_ref"]),
                        "reference_ids": split["public_ref"],
                        "run": lambda problem=problem, public_ref_matrix=public_ref_matrix, beacon_matrix=beacon_matrix, runtime_config=runtime_config, seeds=seeds: _with_config(
                            run_saleem_attack(problem, public_ref_matrix, beacon_matrix, args),
                            config=runtime_config,
                            optimizer_seed=seeds["optimizer_seed"],
                        ),
                    },
                    {
                        "method": "matched_template_n100",
                        "reference_population": "eur",
                        "template_design": "matched",
                        "panel_size": len(panel_ids),
                        "reference_ids": panel_ids,
                        "run": lambda problem=problem, template_matrix=template_matrix, runtime_config=runtime_config, seeds=seeds: _with_config(
                            reference_haplotype_reconstruction(problem, template_matrix, metadata=runtime_config),
                            config=runtime_config,
                            optimizer_seed=seeds["optimizer_seed"],
                        ),
                    },
                ]
                run_specs = [spec for spec in run_specs if spec["method"] in run_methods]

                for spec in run_specs:
                    result, wall_time_sec = run_timed(spec["run"])
                    result = attach_reference_panel_mode(result, shared=True)
                    runtime_fields = runtime_metadata_fields(result, wall_time_sec)
                    source_overlap_count = len(set(split["beacon"]).intersection(spec["reference_ids"]))
                    evaluation = evaluate_reconstruction(
                        result,
                        truth,
                        target_counts=target_counts,
                        reference_pairwise=reference_pairwise,
                        source_sets={"beacon": split["beacon"], "template_reference": spec["reference_ids"]},
                        accepted_manifest_ids=None,
                    )
                    rows.append(
                        _record(
                            snp_count=args.snp_count,
                            snp_set_replicate=snp_set_replicate,
                            beacon_size=beacon_size,
                            replicate_id=replicate_id,
                            method=spec["method"],
                            reference_population=spec["reference_population"],
                            template_design=spec["template_design"],
                            panel_size=spec["panel_size"],
                            selected_keys=selected_keys,
                            result=result,
                            evaluation=evaluation,
                            source_overlap_count=source_overlap_count,
                            runtime_fields=runtime_fields,
                        )
                    )
                print(f"completed {split_id}")

    aggregate = _aggregate(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(args.out_dir / "saleem_random_m1000_anchor_metrics.csv", rows)
    _write_csv(args.out_dir / "saleem_random_m1000_anchor_aggregate.csv", aggregate)
    _write_json(
        args.out_dir / "saleem_random_m1000_anchor_manifest.json",
        {"manifest_id": MANIFEST_ID, "config": config, "snp_sets": snp_set_manifest},
    )
    _write_report(args.out_dir / "saleem_random_m1000_anchor_report.md", aggregate_rows=aggregate, config=config)
    print(f"wrote {args.out_dir / 'saleem_random_m1000_anchor_report.md'}")


if __name__ == "__main__":
    main()
