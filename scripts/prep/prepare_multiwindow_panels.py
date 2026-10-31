"""Prepare an expanded local multi-window panel from existing candidate pfiles."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
PLINK_DIR = DATA_DIR / "plink"
DEFAULT_OUT_DIR = DATA_DIR / "genome_biology_expanded" / "multiwindow"
WRONG_POP = "AFR"


@dataclass(frozen=True)
class WindowSpec:
    chrom: str
    window_index: int
    start: int
    end: int
    span: int
    variant_ids: tuple[str, ...]
    window_size: int

    @property
    def window_id(self) -> str:
        return f"chr{self.chrom}_{self.start}_{self.end}_dense{self.window_size}_w{self.window_index}"


def _candidate_prefix(chrom: str) -> Path:
    return PLINK_DIR / f"1000g_hg38_eur_chr{chrom}_maf01_25_candidates"


def _read_candidate_variants(chrom: str, *, window_size: int) -> list[tuple[str, int]]:
    path = _candidate_prefix(chrom).with_suffix(".pvar")
    rows: list[tuple[str, int]] = []
    seen_positions: set[int] = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                parts = line.rstrip("\n").split()
            row_chrom, pos, variant_id = parts[0], int(parts[1]), parts[2]
            normalized = row_chrom[3:] if row_chrom.startswith("chr") else row_chrom
            if normalized != chrom or pos in seen_positions:
                continue
            seen_positions.add(pos)
            rows.append((variant_id, pos))
    rows.sort(key=lambda item: item[1])
    if len(rows) < window_size:
        raise ValueError(f"chr{chrom} has only {len(rows)} unique-position variants.")
    return rows


def _select_dense_windows(
    chrom: str,
    *,
    window_size: int,
    windows_per_chrom: int,
    min_gap_bp: int,
) -> list[WindowSpec]:
    variants = _read_candidate_variants(chrom, window_size=window_size)
    candidates: list[WindowSpec] = []
    for start_index in range(0, len(variants) - window_size + 1):
        window = variants[start_index : start_index + window_size]
        start = window[0][1]
        end = window[-1][1]
        candidates.append(
            WindowSpec(
                chrom=chrom,
                window_index=-1,
                start=start,
                end=end,
                span=end - start,
                variant_ids=tuple(variant_id for variant_id, _pos in window),
                window_size=window_size,
            )
        )
    selected: list[WindowSpec] = []
    for candidate in sorted(candidates, key=lambda item: (item.span, item.start)):
        if all(abs(candidate.start - existing.start) >= min_gap_bp for existing in selected):
            selected.append(
                WindowSpec(
                    chrom=chrom,
                    window_index=len(selected),
                    start=candidate.start,
                    end=candidate.end,
                    span=candidate.span,
                    variant_ids=candidate.variant_ids,
                    window_size=window_size,
                )
            )
        if len(selected) == windows_per_chrom:
            break
    if len(selected) != windows_per_chrom:
        raise ValueError(f"Could not select {windows_per_chrom} separated windows for chr{chrom}.")
    return selected


def _write_snplist(window: WindowSpec, *, out_dir: Path) -> Path:
    path = out_dir / "snplists" / f"{window.window_id}.snplist"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(window.variant_ids) + "\n")
    return path


def _resolve_tool(explicit: str | None, name: str) -> Path:
    if explicit:
        path = Path(explicit)
        if not path.exists():
            raise FileNotFoundError(f"{name} not found: {path}")
        return path
    found = shutil.which(name)
    if found is None:
        raise FileNotFoundError(f"{name} was not found on PATH. Pass --{name}.")
    return Path(found)


def _run_plink_extract(window: WindowSpec, snplist: Path, *, out_dir: Path, plink2: Path) -> None:
    out_prefix = out_dir / "plink" / window.window_id
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(plink2),
        "--pfile",
        str(_candidate_prefix(window.chrom)),
        "--extract",
        str(snplist),
        "--export",
        "A",
        "--out",
        str(out_prefix),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)


def _write_wrong_population_keep_file(*, out_dir: Path, source_psam: Path) -> Path:
    path = out_dir / "wrong_population" / f"{WRONG_POP.lower()}_samples.psam"
    path.parent.mkdir(parents=True, exist_ok=True)
    with source_psam.open() as source, path.open("w") as out:
        header = source.readline()
        out.write("#IID\n")
        columns = header.strip().lstrip("#").split("\t")
        iid_index = columns.index("IID")
        super_pop_index = columns.index("SuperPop")
        for line in source:
            parts = line.rstrip("\n").split("\t")
            if parts[super_pop_index] == WRONG_POP:
                out.write(parts[iid_index] + "\n")
    return path


def _run_wrong_population_export(
    window: WindowSpec,
    snplist: Path,
    keep_file: Path,
    *,
    out_dir: Path,
    plink2: Path,
    source_pfile: Path,
) -> None:
    out_prefix = out_dir / "wrong_population" / window.window_id
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(plink2),
        "--pfile",
        str(source_pfile),
        "vzs",
        "--allow-extra-chr",
        "--keep",
        str(keep_file),
        "--extract",
        str(snplist),
        "--export",
        "A",
        "--out",
        str(out_prefix),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chroms", default="1,20,22")
    parser.add_argument("--window-size", type=int, default=500)
    parser.add_argument("--windows-per-chrom", type=int, default=5)
    parser.add_argument("--min-gap-bp", type=int, default=1_000_000)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--plink2", default=None, help="Path to plink2. Defaults to plink2 on PATH.")
    parser.add_argument("--pfile", type=Path, required=True, help="1000 Genomes pfile prefix for the full callset.")
    parser.add_argument("--psam", type=Path, default=None, help="Sample file. Defaults to <pfile>.psam.")
    args = parser.parse_args()
    plink2 = _resolve_tool(args.plink2, "plink2")
    source_pfile = args.pfile
    source_psam = args.psam if args.psam is not None else source_pfile.with_suffix(".psam")
    if not source_psam.exists():
        raise FileNotFoundError(f"Sample file not found: {source_psam}")

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    chroms = [chrom.strip() for chrom in args.chroms.split(",") if chrom.strip()]
    wrong_pop_keep = _write_wrong_population_keep_file(out_dir=out_dir, source_psam=source_psam)
    all_windows: list[WindowSpec] = []
    for chrom in chroms:
        windows = _select_dense_windows(
            chrom,
            window_size=args.window_size,
            windows_per_chrom=args.windows_per_chrom,
            min_gap_bp=args.min_gap_bp,
        )
        for window in windows:
            snplist = _write_snplist(window, out_dir=out_dir)
            _run_plink_extract(window, snplist, out_dir=out_dir, plink2=plink2)
            _run_wrong_population_export(
                window,
                snplist,
                wrong_pop_keep,
                out_dir=out_dir,
                plink2=plink2,
                source_pfile=source_pfile,
            )
            all_windows.append(window)
            print(f"prepared {window.window_id} span={window.span}")

    manifest_path = out_dir / "multiwindow_manifest.csv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "window_id",
                "chrom",
                "window_index",
                "start",
                "end",
                "span",
                "variant_count",
                "raw_path",
                "snplist_path",
            ],
        )
        writer.writeheader()
        for window in all_windows:
            writer.writerow(
                {
                    "window_id": window.window_id,
                    "chrom": window.chrom,
                    "window_index": window.window_index,
                    "start": window.start,
                    "end": window.end,
                    "span": window.span,
                    "variant_count": len(window.variant_ids),
                    "raw_path": str(out_dir / "plink" / f"{window.window_id}.raw"),
                    "snplist_path": str(out_dir / "snplists" / f"{window.window_id}.snplist"),
                }
            )
    print(f"wrote {manifest_path}")


if __name__ == "__main__":
    main()
