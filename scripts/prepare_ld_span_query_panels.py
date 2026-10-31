"""Prepare fixed-size SNP panels that vary LD, density, and genomic span.

The main analysis uses six compact dense windows. This constructor adds
further regions on the same candidate chromosomes and, at each region, builds
six panels with the same SNP count:

- dense: consecutive candidates
- ld_r2_080 / ld_r2_050 / ld_r2_010: weaker residual linkage
- sparse: minimum gap of 50 kb
- long_range: even spacing across 10 Mb

LD pruning uses carrier genotypes from the EUR candidate panel, matching the
existing local query-set constructor. Selected panels are disjoint from the
six main dense windows.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
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

from beacon_reconstruction.query_regimes import (
    REGIME_SPECS,
    RegimeSpec,
    choose_anchor_indices,
    correlation_summaries,
    eligible_indices,
    greedy_ld_prune,
    panel_geometry,
    select_consecutive,
    select_even_span,
    select_target_gap,
)
from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers


def _resolve_plink2(explicit: str | None) -> Path:
    if explicit:
        path = Path(explicit)
        if not path.exists():
            raise FileNotFoundError(f"plink2 not found: {path}")
        return path
    found = shutil.which("plink2")
    if found is None:
        raise FileNotFoundError(
            "plink2 was not found on PATH. Pass --plink2. "
            "Prepared panels are already in data/ld_span_queries and do not need this step."
        )
    return Path(found)
PLINK_DIR = ROOT / "data" / "plink"
MAIN_MANIFEST = ROOT / "data" / "multiwindow" / "multiwindow_manifest.csv"
DEFAULT_OUT_DIR = ROOT / "data" / "ld_span_queries"
CHROMS = ("1", "20", "22")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _candidate_prefix(chrom: str) -> Path:
    return PLINK_DIR / f"1000g_hg38_eur_chr{chrom}_maf01_25_candidates"


def _read_candidate_variants(chrom: str) -> list[tuple[str, int]]:
    path = _candidate_prefix(chrom).with_suffix(".pvar")
    rows: list[tuple[str, int]] = []
    seen_positions: set[int] = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom_field, pos_field, variant_id, _rest = line.split("\t", 3)
            normalized = chrom_field[3:] if chrom_field.startswith("chr") else chrom_field
            pos = int(pos_field)
            if normalized != chrom or pos in seen_positions:
                continue
            seen_positions.add(pos)
            rows.append((variant_id, pos))
    rows.sort(key=lambda item: item[1])
    if not rows:
        raise ValueError(f"No candidate variants for chr{chrom}.")
    return rows


def _excluded_intervals(chrom: str, *, flank_bp: int) -> list[tuple[int, int]]:
    if not MAIN_MANIFEST.exists():
        return []
    intervals: list[tuple[int, int]] = []
    for row in _read_csv(MAIN_MANIFEST):
        if row["chrom"] != chrom:
            continue
        intervals.append((int(row["start"]) - flank_bp, int(row["end"]) + flank_bp))
    return intervals


def _write_snplist(path: Path, variant_ids: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(variant_ids) + "\n")


def _export_raw(chrom: str, snplist: Path, out_prefix: Path, *, force: bool, plink2: Path) -> Path:
    raw_path = Path(str(out_prefix) + ".raw")
    if raw_path.exists() and not force:
        return raw_path
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(plink2),
        "--pfile",
        str(_candidate_prefix(chrom)),
        "--extract",
        str(snplist),
        "--export",
        "A",
        "--out",
        str(out_prefix),
    ]
    completed = subprocess.run(cmd, cwd=ROOT, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"plink2 failed for {out_prefix} (exit {completed.returncode}).\n"
            f"{completed.stdout[-2000:]}\n{completed.stderr[-2000:]}"
        )
    if not raw_path.exists():
        raise FileNotFoundError(f"plink2 did not write {raw_path}.")
    return raw_path


def _remove_plink_outputs(prefix: Path) -> None:
    for path in prefix.parent.glob(prefix.name + ".*"):
        path.unlink()


def _select_regime_indices(
    spec: RegimeSpec,
    *,
    positions: np.ndarray,
    anchor: int,
    snp_count: int,
    excluded: list[tuple[int, int]],
    ld_pool_span_bp: int,
    dosage_for_pool: np.ndarray | None,
    pool_indices: np.ndarray | None,
) -> np.ndarray:
    if spec.regime == "dense":
        return select_consecutive(anchor, snp_count, len(positions))
    if spec.target_gap_bp is not None:
        pool = eligible_indices(
            positions,
            start_index=anchor,
            end_pos=None,
            excluded=excluded,
        )
        return select_target_gap(positions, pool, count=snp_count, target_gap_bp=spec.target_gap_bp)
    if spec.target_span_bp is not None:
        pool = eligible_indices(
            positions,
            start_index=anchor,
            end_pos=int(positions[anchor]) + spec.target_span_bp,
            excluded=excluded,
        )
        return select_even_span(positions, pool, count=snp_count, target_span_bp=spec.target_span_bp)
    if spec.r2_threshold is None:
        raise ValueError(f"Regime {spec.regime} has no selection rule.")
    if dosage_for_pool is None or pool_indices is None:
        raise ValueError(f"Regime {spec.regime} requires a genotype pool.")
    local = greedy_ld_prune(dosage_for_pool, r2_threshold=spec.r2_threshold, count=snp_count)
    if len(local) != snp_count:
        raise ValueError(
            f"Regime {spec.regime} kept {len(local)} variants at r^2<={spec.r2_threshold}; needed {snp_count}."
        )
    return pool_indices[local]


def _load_pool_dosage(raw_path: Path, variant_ids: list[str]) -> np.ndarray:
    _sample_ids, loaded_ids, _positions, _targets, carriers = _load_raw_as_target_carriers(raw_path)
    if tuple(loaded_ids) != tuple(variant_ids):
        raise ValueError(f"Pool variant order mismatch for {raw_path}.")
    return carriers


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chroms", default=",".join(CHROMS))
    parser.add_argument("--regions-per-chrom", type=int, default=2)
    parser.add_argument("--snp-count", type=int, default=200)
    parser.add_argument("--min-gap-bp", type=int, default=12_000_000)
    parser.add_argument("--exclude-flank-bp", type=int, default=500_000)
    parser.add_argument("--ld-pool-span-bp", type=int, default=30_000_000)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--plink2", default=None, help="Path to plink2. Defaults to plink2 on PATH.")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    plink2 = _resolve_plink2(args.plink2)

    chroms = [value.strip() for value in args.chroms.split(",") if value.strip()]
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []

    for chrom in chroms:
        variants = _read_candidate_variants(chrom)
        variant_ids = [variant_id for variant_id, _pos in variants]
        positions = np.asarray([pos for _variant_id, pos in variants], dtype=np.int64)
        excluded = _excluded_intervals(chrom, flank_bp=args.exclude_flank_bp)
        sparse_gap = next(spec.target_gap_bp for spec in REGIME_SPECS if spec.regime == "sparse")
        long_span = next(spec.target_span_bp for spec in REGIME_SPECS if spec.regime == "long_range")
        assert sparse_gap is not None and long_span is not None
        anchors = choose_anchor_indices(
            positions,
            anchors=args.regions_per_chrom,
            snp_count=args.snp_count,
            min_gap_bp=args.min_gap_bp,
            excluded=excluded,
            long_range_span_bp=long_span,
            sparse_gap_bp=sparse_gap,
            ld_pool_span_bp=args.ld_pool_span_bp,
            min_following_variants=args.snp_count * 20,
        )
        print(f"chr{chrom}: {len(variants)} candidates, anchors={anchors}")

        for region_index, anchor in enumerate(anchors):
            pool_end = int(positions[anchor]) + args.ld_pool_span_bp
            pool_indices = eligible_indices(
                positions,
                start_index=anchor,
                end_pos=pool_end,
                excluded=excluded,
            )
            pool_id = f"chr{chrom}_anchor{region_index}_ldpool"
            pool_snplist = out_dir / "pools" / f"{pool_id}.snplist"
            pool_ids = [variant_ids[int(index)] for index in pool_indices]
            _write_snplist(pool_snplist, pool_ids)
            pool_prefix = out_dir / "pools" / pool_id
            pool_raw = _export_raw(chrom, pool_snplist, pool_prefix, force=args.force, plink2=plink2)
            pool_dosage = _load_pool_dosage(pool_raw, pool_ids)

            for spec in REGIME_SPECS:
                selected = _select_regime_indices(
                    spec,
                    positions=positions,
                    anchor=anchor,
                    snp_count=args.snp_count,
                    excluded=excluded,
                    ld_pool_span_bp=args.ld_pool_span_bp,
                    dosage_for_pool=pool_dosage,
                    pool_indices=pool_indices,
                )
                selected_positions = positions[selected]
                geometry = panel_geometry(selected_positions)
                query_set_id = (
                    f"chr{chrom}_{int(geometry['start'])}_{int(geometry['end'])}_{spec.regime}_r{region_index}"
                )
                snplist_path = out_dir / "snplists" / f"{query_set_id}.snplist"
                selected_ids = [variant_ids[int(index)] for index in selected]
                _write_snplist(snplist_path, selected_ids)
                raw_path = _export_raw(
                    chrom,
                    snplist_path,
                    out_dir / "plink" / query_set_id,
                    force=args.force,
                    plink2=plink2,
                )
                _samples, loaded_ids, _pos, _targets, carriers = _load_raw_as_target_carriers(raw_path)
                if tuple(loaded_ids) != tuple(selected_ids):
                    raise ValueError(f"Exported variant order mismatch for {query_set_id}.")
                summaries = correlation_summaries(carriers)
                rows.append(
                    {
                        "query_set_id": query_set_id,
                        "regime": spec.regime,
                        "description": spec.description,
                        "chrom": chrom,
                        "region_index": region_index,
                        "anchor_position": int(positions[anchor]),
                        "region_role": "additional",
                        "snp_count_requested": args.snp_count,
                        "variant_count": int(len(selected)),
                        "start": int(geometry["start"]),
                        "end": int(geometry["end"]),
                        "span_bp": int(geometry["span_bp"]),
                        "median_inter_snp_bp": geometry["median_inter_snp_bp"],
                        "snps_per_mb": geometry["snps_per_mb"],
                        "mean_abs_corr": summaries["mean_abs_corr"],
                        "max_pairwise_r2": summaries["max_pairwise_r2"],
                        "r2_threshold": "" if spec.r2_threshold is None else spec.r2_threshold,
                        "target_gap_bp": "" if spec.target_gap_bp is None else spec.target_gap_bp,
                        "target_span_bp": "" if spec.target_span_bp is None else spec.target_span_bp,
                        "raw_path": str(raw_path),
                        "snplist_path": str(snplist_path),
                    }
                )
                print(
                    f"prepared {query_set_id} span_mb={geometry['span_bp'] / 1e6:.3f} "
                    f"median_gap={geometry['median_inter_snp_bp']:.0f} "
                    f"|r|={summaries['mean_abs_corr']:.3f}"
                )
            _remove_plink_outputs(pool_prefix)

    manifest_path = out_dir / "ld_span_query_manifest.csv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {manifest_path} ({len(rows)} query sets)")


if __name__ == "__main__":
    main()
