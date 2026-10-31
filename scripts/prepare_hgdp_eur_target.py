"""Slice gnomAD HGDP European genotypes onto the existing hg38 windows.

The beacon cohort is the HGDP European subset of the gnomAD v3.1.2 HGDP+1KG
callset. 1000 Genomes samples in that callset are not used. Alleles are coded
with the same target allele as the local 1000 Genomes window.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
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

from scripts.lib.multiwindow_data import load_raw_carriers as _load_raw_as_target_carriers


DEFAULT_MULTIWINDOW_DIR = ROOT / "data" / "multiwindow"
DEFAULT_OUT_DIR = ROOT / "data" / "external" / "hgdp_1kg"
MIN_SAMPLE_CALL_RATE = 0.95
METADATA_NAME = "gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz"
VCF_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/"
    "gnomad.genomes.v3.1.2.hgdp_tgp.chr{chrom}.vcf.bgz"
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def load_hgdp_eur_metadata(metadata_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with gzip.open(metadata_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: position for position, name in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            meta = json.loads(parts[index["hgdp_tgp_meta"]])
            if meta.get("project") != "HGDP" or meta.get("genetic_region") != "EUR":
                continue
            rows.append(
                {
                    "sample_id": parts[0],
                    "population": str(meta.get("population") or ""),
                    "genetic_region": "EUR",
                    "project": "HGDP",
                    "gnomad_high_quality": parts[index["gnomad_high_quality"]],
                    "high_quality": parts[index["high_quality"]],
                }
            )
    rows.sort(key=lambda row: row["sample_id"])
    if not rows:
        raise ValueError(f"No HGDP European samples in {metadata_path}.")
    return rows


def _allele_dosage(genotype: str, allele_index: int) -> int | None:
    if genotype in {"", ".", "./.", ".|."}:
        return None
    separator = "|" if "|" in genotype else "/"
    called = genotype.split(separator)
    if any(part == "." for part in called):
        return None
    try:
        return sum(int(part) == allele_index for part in called)
    except ValueError:
        return None


def _choose_record(
    records: list[dict[str, Any]],
    *,
    ref: str,
    alt: str,
) -> tuple[dict[str, Any] | None, str]:
    exact = [record for record in records if record["ref"] == ref and record["alt"] == alt]
    if exact:
        return exact[0], "exact"
    flipped = [record for record in records if record["ref"] == alt and record["alt"] == ref]
    if flipped:
        return flipped[0], "flipped"
    containing = [
        record
        for record in records
        if record["ref"] == ref and alt in record["alt"].split(",")
    ]
    if containing:
        return containing[0], "multiallelic"
    return None, "missing"


def retain_called_subset(
    missing: np.ndarray,
    *,
    min_call_rate: float,
) -> np.ndarray:
    """Return samples whose call rate is at least ``min_call_rate``.

    ``missing`` covers allele-matched sites only, with one row per sample.
    """

    if missing.ndim != 2 or missing.shape[1] == 0:
        raise ValueError("missing must be a samples-by-variants matrix with at least one variant.")
    call_rate = 1.0 - missing.mean(axis=1)
    kept = call_rate >= min_call_rate
    if not np.any(kept):
        raise ValueError(f"No sample reaches call rate {min_call_rate}.")
    return kept


def harmonize_genotypes(
    *,
    variant_ids: tuple[str, ...],
    target_alleles: dict[str, str],
    records_by_position: dict[int, list[dict[str, Any]]],
    sample_count: int,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    carriers = np.zeros((sample_count, len(variant_ids)), dtype=np.int8)
    status_rows: list[dict[str, Any]] = []
    for column, variant_id in enumerate(variant_ids):
        _chrom, position_text, ref, alt = variant_id.split(":")
        position = int(position_text)
        record, match = _choose_record(records_by_position.get(position, []), ref=ref, alt=alt)
        target_allele = target_alleles[variant_id]
        row = {
            "variant_id": variant_id,
            "match": match,
            "target_allele": target_allele,
            "kept": False,
            "missing_genotypes": "",
        }
        if record is None:
            status_rows.append(row)
            continue
        alleles = [record["ref"], *record["alt"].split(",")]
        if target_allele not in alleles:
            row["match"] = "target_allele_absent"
            status_rows.append(row)
            continue
        allele_index = alleles.index(target_allele)
        column_carriers = np.zeros(sample_count, dtype=np.int8)
        missing = 0
        for sample_index, genotype in enumerate(record["genotypes"]):
            dosage = _allele_dosage(genotype, allele_index)
            if dosage is None:
                missing += 1
                continue
            column_carriers[sample_index] = int(dosage > 0)
        row["missing_genotypes"] = missing
        if missing:
            row["match"] = f"{match}_missing"
            status_rows.append(row)
            continue
        carriers[:, column] = column_carriers
        row["kept"] = True
        status_rows.append(row)
    return carriers, status_rows


def build_call_matrix(
    *,
    variant_ids: tuple[str, ...],
    target_alleles: dict[str, str],
    records_by_position: dict[int, list[dict[str, Any]]],
    sample_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Return carriers, a missing mask, an allele-matched mask, and status rows."""

    carriers = np.zeros((sample_count, len(variant_ids)), dtype=np.int8)
    missing = np.ones((sample_count, len(variant_ids)), dtype=bool)
    matched = np.zeros(len(variant_ids), dtype=bool)
    status_rows: list[dict[str, Any]] = []
    for column, variant_id in enumerate(variant_ids):
        _chrom, position_text, ref, alt = variant_id.split(":")
        record, match = _choose_record(records_by_position.get(int(position_text), []), ref=ref, alt=alt)
        target_allele = target_alleles[variant_id]
        row = {
            "variant_id": variant_id,
            "match": match,
            "target_allele": target_allele,
            "kept": False,
            "missing_genotypes": "",
        }
        if record is None:
            status_rows.append(row)
            continue
        alleles = [record["ref"], *record["alt"].split(",")]
        if target_allele not in alleles:
            row["match"] = "target_allele_absent"
            status_rows.append(row)
            continue
        allele_index = alleles.index(target_allele)
        missing_count = 0
        for sample_index, genotype in enumerate(record["genotypes"]):
            dosage = _allele_dosage(genotype, allele_index)
            if dosage is None:
                missing_count += 1
                continue
            missing[sample_index, column] = False
            carriers[sample_index, column] = int(dosage > 0)
        matched[column] = True
        row["missing_genotypes"] = missing_count
        status_rows.append(row)
    return carriers, missing, matched, status_rows


def _query_region(
    *,
    chrom: str,
    start: int,
    end: int,
    sample_file: Path,
    destination: Path,
) -> None:
    if destination.exists() and destination.stat().st_size > 0:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    url = VCF_URL.format(chrom=chrom)
    command = [
        "bcftools",
        "query",
        "-r",
        f"chr{chrom}:{start}-{end}",
        "-S",
        str(sample_file),
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
        url,
    ]
    last_error = ""
    for attempt in range(1, 4):
        with destination.open("w") as handle:
            completed = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True)
        if completed.returncode == 0 and destination.stat().st_size > 0:
            return
        last_error = completed.stderr.strip()
        if completed.returncode == 0 and "error" not in last_error.lower() and "fail" not in last_error.lower():
            destination.write_text("# no_variants\n")
            print(f"no gnomAD variants in chr{chrom}:{start}-{end}", flush=True)
            return
        destination.write_text("")
        print(f"retry {attempt} for chr{chrom}:{start}-{end}: {last_error[-400:]}", flush=True)
    raise RuntimeError(f"Empty genotype slice for chr{chrom}:{start}-{end}. {last_error}")


def _parse_query(path: Path, sample_count: int) -> dict[int, list[dict[str, Any]]]:
    records: dict[int, list[dict[str, Any]]] = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            genotypes = parts[4:]
            if len(genotypes) != sample_count:
                raise ValueError(
                    f"{path} has {len(genotypes)} genotypes at {parts[0]}:{parts[1]}, expected {sample_count}."
                )
            position = int(parts[1])
            records.setdefault(position, []).append(
                {"ref": parts[2], "alt": parts[3], "genotypes": genotypes}
            )
    return records


def _vcf_sample_order(chrom: str, cache_path: Path) -> list[str]:
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return [line.strip() for line in cache_path.read_text().splitlines() if line.strip()]
    command = ["bcftools", "query", "-l", VCF_URL.format(chrom=chrom)]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    sample_ids = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if not sample_ids:
        raise ValueError(f"No samples in the gnomAD header for chromosome {chrom}.")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text("".join(f"{sample_id}\n" for sample_id in sample_ids))
    return sample_ids


def prepare_windows(*, multiwindow_dir: Path, out_dir: Path, max_windows: int) -> None:
    metadata_rows = load_hgdp_eur_metadata(out_dir / METADATA_NAME)
    by_id = {row["sample_id"]: row for row in metadata_rows}
    windows = _read_csv(multiwindow_dir / "multiwindow_manifest.csv")[: max(0, max_windows)]
    vcf_order = _vcf_sample_order(windows[0]["chrom"], out_dir / "hgdp_tgp_vcf_sample_order.txt")
    metadata_rows = [by_id[sample_id] for sample_id in vcf_order if sample_id in by_id]
    missing = sorted(set(by_id).difference(row["sample_id"] for row in metadata_rows))
    if missing:
        raise ValueError(f"{len(missing)} HGDP European samples are absent from the VCF header.")
    sample_ids = [row["sample_id"] for row in metadata_rows]
    _write_csv(out_dir / "hgdp_eur_samples.csv", metadata_rows)
    sample_file = out_dir / "hgdp_eur_samples.txt"
    sample_file.write_text("".join(f"{sample_id}\n" for sample_id in sample_ids))
    window_dir = out_dir / "windows"
    parsed: list[dict[str, Any]] = []
    for window in windows:
        window_id = window["window_id"]
        _sample_ids_1kg, variant_ids, positions, target_alleles, _matrix = _load_raw_as_target_carriers(
            Path(window["raw_path"])
        )
        query_path = window_dir / f"{window_id}.genotypes.tsv"
        print(f"query {window_id}", flush=True)
        _query_region(
            chrom=window["chrom"],
            start=int(window["start"]),
            end=int(window["end"]),
            sample_file=sample_file,
            destination=query_path,
        )
        records = _parse_query(query_path, len(sample_ids))
        carriers, missing_mask, matched, status_rows = build_call_matrix(
            variant_ids=variant_ids,
            target_alleles=target_alleles,
            records_by_position=records,
            sample_count=len(sample_ids),
        )
        parsed.append(
            {
                "window_id": window_id,
                "variant_ids": variant_ids,
                "positions": positions,
                "carriers": carriers,
                "missing": missing_mask,
                "matched": matched,
                "status_rows": status_rows,
            }
        )

    matched_blocks = [item["missing"][:, item["matched"]] for item in parsed if np.any(item["matched"])]
    if not matched_blocks:
        raise ValueError("No HGDP sites matched the 1000 Genomes target alleles.")
    retained = retain_called_subset(np.concatenate(matched_blocks, axis=1), min_call_rate=MIN_SAMPLE_CALL_RATE)
    if int(np.sum(retained)) < 100:
        raise ValueError(
            f"Only {int(np.sum(retained))} HGDP samples reach call rate {MIN_SAMPLE_CALL_RATE}."
        )
    joint_missing = np.concatenate(matched_blocks, axis=1)
    call_rate = 1.0 - joint_missing.mean(axis=1)
    for row, keep, rate in zip(metadata_rows, retained, call_rate):
        row["retained"] = bool(keep)
        row["call_rate"] = f"{float(rate):.6f}"
    _write_csv(out_dir / "hgdp_eur_samples.csv", metadata_rows)
    retained_ids = [sample_id for sample_id, keep in zip(sample_ids, retained) if keep]
    print(
        f"retained {len(retained_ids)}/{len(sample_ids)} HGDP EUR samples at call rate >= {MIN_SAMPLE_CALL_RATE}",
        flush=True,
    )

    summary: list[dict[str, Any]] = []
    for item in parsed:
        window_id = item["window_id"]
        status_rows = item["status_rows"]
        complete = item["matched"] & ~item["missing"][retained].any(axis=0)
        kept_index = np.flatnonzero(complete)
        for index, row in enumerate(status_rows):
            row["kept"] = bool(complete[index])
            if item["matched"][index] and not complete[index]:
                row["match"] = f"{row['match']}_missing"
        kept_carriers = item["carriers"][np.ix_(retained, kept_index)]
        np.savez_compressed(
            window_dir / f"{window_id}.npz",
            sample_ids=np.asarray(retained_ids),
            populations=np.asarray([row["population"] for row, keep in zip(metadata_rows, retained) if keep]),
            variant_ids=np.asarray([item["variant_ids"][index] for index in kept_index]),
            positions=np.asarray([item["positions"][index] for index in kept_index], dtype=np.int64),
            carriers=kept_carriers,
        )
        _write_csv(window_dir / f"{window_id}.harmonization.csv", status_rows)
        summary.append(
            {
                "window_id": window_id,
                "target_variants": len(item["variant_ids"]),
                "kept_variants": int(kept_index.size),
                "exact": sum(row["match"] == "exact" and row["kept"] for row in status_rows),
                "flipped": sum(row["match"] == "flipped" and row["kept"] for row in status_rows),
                "multiallelic": sum(row["match"] == "multiallelic" and row["kept"] for row in status_rows),
                "missing_site": sum(row["match"] == "missing" for row in status_rows),
                "dropped_missing_genotype": sum(str(row["match"]).endswith("_missing") for row in status_rows),
                "hgdp_samples": len(retained_ids),
                "min_sample_call_rate": MIN_SAMPLE_CALL_RATE,
            }
        )
        print(
            f"kept {window_id}: {int(kept_index.size)}/{len(item['variant_ids'])} variants, "
            f"{len(retained_ids)} HGDP EUR samples",
            flush=True,
        )
    _write_csv(out_dir / "hgdp_eur_window_overlap.csv", summary)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--multiwindow-dir", type=Path, default=DEFAULT_MULTIWINDOW_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--max-windows", type=int, default=6)
    args = parser.parse_args()
    prepare_windows(multiwindow_dir=args.multiwindow_dir, out_dir=args.out_dir, max_windows=args.max_windows)


if __name__ == "__main__":
    main()
