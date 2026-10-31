"""Launch explicit target x wrong-ancestry reference-matching runs.

This wrapper keeps ``run_multiancestry_reference_matching.py`` unchanged.  The
runner currently chooses one internal wrong-reference label per target, so this
script creates a small per-comparison panel manifest that aliases that internal
label to the requested wrong ancestry, runs the existing CLI, and then rewrites
the output labels to the real requested ancestry.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PANEL_DIR = ROOT / "data" / "multiwindow" / "ancestry_panels"
DEFAULT_OUT_DIR = ROOT / "results" / "experiments" / "multiancestry_reference_matching_all_wrong"
DEFAULT_SHARED_OUT_DIR = (
    ROOT / "results" / "experiments" / "multiancestry_reference_matching_all_wrong_shared_panel"
)
RUNNER = ROOT / "scripts" / "run_multiancestry_reference_matching.py"
METRICS = "multiancestry_reference_matching_metrics.csv"
COMBINED_METRICS = "all_wrong_combined_metrics.csv"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.run_multiancestry_reference_matching import saleem_methods_for_fractions
from scripts.ld_af_helpers import add_saleem_arguments, merge_saleem_rows
AGGREGATE = "multiancestry_reference_matching_aggregate.csv"
MANIFEST = "multiancestry_reference_matching_manifest.json"
REPORT = "multiancestry_reference_matching_report.md"


def _split_csv(value: str) -> list[str]:
    return [item.strip().upper() for item in value.split(",") if item.strip()]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], *, fieldnames: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}.")
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _runner_wrong_label(target_ancestry: str) -> str:
    return "AFR" if target_ancestry == "EUR" else "EUR"


def _comparison_pairs(target_ancestries: list[str], wrong_ancestries: list[str]) -> list[tuple[str, str]]:
    return [(target, wrong) for target in target_ancestries for wrong in wrong_ancestries if target != wrong]


def _write_aliased_manifest(
    *,
    source_rows: list[dict[str, str]],
    panel_dir: Path,
    target_ancestry: str,
    wrong_ancestry: str,
) -> tuple[Path, str]:
    runner_wrong = _runner_wrong_label(target_ancestry)
    selected: list[dict[str, str]] = []
    for row in source_rows:
        ancestry = row["ancestry"].upper()
        if ancestry == target_ancestry:
            selected.append(dict(row, ancestry=target_ancestry))
        elif ancestry == wrong_ancestry:
            selected.append(dict(row, ancestry=runner_wrong))

    if not selected:
        raise ValueError(f"No panel rows found for {target_ancestry} vs {wrong_ancestry}.")

    required = {target_ancestry, runner_wrong}
    observed = {row["ancestry"] for row in selected}
    missing = required - observed
    if missing:
        raise ValueError(
            f"Panel manifest for {target_ancestry} vs {wrong_ancestry} is missing aliased rows: {sorted(missing)}"
        )

    manifest_path = panel_dir / "ancestry_panel_manifest.csv"
    _write_csv(manifest_path, selected, fieldnames=list(source_rows[0].keys()))
    return manifest_path, runner_wrong


def _replace_reference_label(value: str, *, target: str, runner_wrong: str, actual_wrong: str) -> str:
    if value == runner_wrong:
        return actual_wrong
    mixed_alias = f"{target}+{runner_wrong}"
    if value == mixed_alias:
        return f"{target}+{actual_wrong}"
    return value


def _normalize_csv(path: Path, *, target: str, runner_wrong: str, actual_wrong: str) -> None:
    rows = _read_csv(path)
    for row in rows:
        if "reference_ancestry" in row:
            row["reference_ancestry"] = _replace_reference_label(
                row["reference_ancestry"],
                target=target,
                runner_wrong=runner_wrong,
                actual_wrong=actual_wrong,
            )
    _write_csv(path, rows, fieldnames=list(rows[0].keys()))


def _normalize_manifest(
    path: Path,
    *,
    target: str,
    runner_wrong: str,
    actual_wrong: str,
    aliased_manifest: Path,
    command: list[str],
) -> None:
    payload = json.loads(path.read_text())
    config = payload.setdefault("config", {})
    config["actual_wrong_ancestry"] = actual_wrong
    config["runner_wrong_label"] = runner_wrong
    config["target_wrong_pair"] = {"target_ancestry": target, "wrong_ancestry": actual_wrong}
    config["aliased_panel_manifest"] = str(aliased_manifest)
    config["launcher_command"] = command
    if "reference_panel_mode" not in config:
        config["reference_panel_mode"] = "shared" if "--shared-reference-panel" in command else "dual_disjoint"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _runner_saleem_cli_args(args: argparse.Namespace, *, pair_metrics: Path | None) -> list[str]:
    cli = [
        "--saleem-corr-epoch",
        str(args.saleem_corr_epoch),
        "--saleem-freq-epoch",
        str(args.saleem_freq_epoch),
        "--saleem-max-iterations",
        str(args.saleem_max_iterations),
        "--saleem-flip-threshold",
        str(args.saleem_flip_threshold),
        "--saleem-min-outer-iterations",
        str(args.saleem_min_outer_iterations),
        "--saleem-learning-rate",
        str(args.saleem_learning_rate),
    ]
    if args.shared_reference_panel:
        cli.append("--shared-reference-panel")
    if args.run_saleem:
        cli.append("--run-saleem")
    else:
        cli.append("--no-run-saleem")
    if args.saleem_only:
        cli.append("--saleem-only")
        if pair_metrics is not None:
            cli.extend(["--merge-metrics", str(pair_metrics)])
    return cli


def _annotate_pair_rows(rows: list[dict[str, str]], *, target: str, wrong: str) -> list[dict[str, Any]]:
    return [
        {
            **row,
            "wrong_ancestry": wrong,
            "target_wrong_pair": f"{target}_vs_{wrong}",
        }
        for row in rows
    ]


def _all_completed_pairs(out_dir: Path) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for path in sorted(out_dir.glob("target_*__wrong_*")):
        if not path.is_dir() or not (path / METRICS).exists():
            continue
        stem = path.name.removeprefix("target_")
        target, wrong = stem.split("__wrong_", maxsplit=1)
        pairs.append((target, wrong))
    return pairs


def _rebuild_combined_metrics(
    *,
    args: argparse.Namespace,
    pairs: list[tuple[str, str]] | None = None,
) -> Path | None:
    pair_list = _all_completed_pairs(args.out_dir) if args.saleem_only else (pairs or [])
    if not pair_list and pairs:
        pair_list = pairs
    combined_rows: list[dict[str, Any]] = []
    for target, wrong in pair_list:
        metrics_path = args.out_dir / f"target_{target}__wrong_{wrong}" / METRICS
        if not metrics_path.exists():
            continue
        combined_rows.extend(_annotate_pair_rows(_read_csv(metrics_path), target=target, wrong=wrong))
    if not combined_rows:
        return None
    if args.saleem_only:
        matched_fractions = sorted({int(value) for value in args.matched_fractions.split(",") if value.strip()})
        saleem_methods = saleem_methods_for_fractions(matched_fractions)
        new_saleem_rows = [
            row for row in combined_rows if row["method"] in saleem_methods
        ]
        combined_rows = merge_saleem_rows(
            new_rows=new_saleem_rows,
            merge_metrics=args.merge_metrics,
            saleem_only=True,
            saleem_methods=saleem_methods,
        )
    combined_path = args.out_dir / COMBINED_METRICS
    _write_csv(combined_path, combined_rows)
    return combined_path


def _normalize_report(path: Path, *, target: str, runner_wrong: str, actual_wrong: str) -> None:
    text = path.read_text()
    text = text.replace(f"| {target} | {runner_wrong} |", f"| {target} | {actual_wrong} |")
    text = text.replace(f"| {target} | {target}+{runner_wrong} |", f"| {target} | {target}+{actual_wrong} |")
    path.write_text(text)


def _run_comparison(
    *,
    args: argparse.Namespace,
    source_rows: list[dict[str, str]],
    target: str,
    wrong: str,
) -> dict[str, Any]:
    pair_name = f"target_{target}__wrong_{wrong}"
    pair_out_dir = args.out_dir / pair_name
    aliased_panel_dir = args.out_dir / "_aliased_panels" / pair_name
    complete_marker = pair_out_dir / METRICS
    if complete_marker.exists() and not args.overwrite:
        return {
            "target_ancestry": target,
            "wrong_ancestry": wrong,
            "status": "skipped_existing",
            "out_dir": str(pair_out_dir),
            "aliased_manifest": str(aliased_panel_dir / "ancestry_panel_manifest.csv"),
        }

    aliased_manifest, runner_wrong = _write_aliased_manifest(
        source_rows=source_rows,
        panel_dir=aliased_panel_dir,
        target_ancestry=target,
        wrong_ancestry=wrong,
    )
    command = [
        sys.executable,
        str(RUNNER),
        "--panel-dir",
        str(aliased_panel_dir),
        "--out-dir",
        str(pair_out_dir),
        "--target-ancestries",
        target,
        "--beacon-sizes",
        args.beacon_sizes,
        "--max-windows",
        str(args.max_windows),
        "--replicates",
        str(args.replicates),
        "--matched-fractions",
        args.matched_fractions,
        "--epochs",
        str(args.epochs),
        "--proposals-per-epoch",
        str(args.proposals_per_epoch),
        "--early-stopping",
        str(args.early_stopping),
        *_runner_saleem_cli_args(
            args,
            pair_metrics=pair_out_dir / METRICS if args.saleem_only else None,
        ),
    ]

    if args.dry_run:
        status = "dry_run"
    else:
        subprocess.run(command, cwd=ROOT, check=True)
        _normalize_csv(pair_out_dir / METRICS, target=target, runner_wrong=runner_wrong, actual_wrong=wrong)
        _normalize_csv(pair_out_dir / AGGREGATE, target=target, runner_wrong=runner_wrong, actual_wrong=wrong)
        _normalize_manifest(
            pair_out_dir / MANIFEST,
            target=target,
            runner_wrong=runner_wrong,
            actual_wrong=wrong,
            aliased_manifest=aliased_manifest,
            command=command,
        )
        _normalize_report(pair_out_dir / REPORT, target=target, runner_wrong=runner_wrong, actual_wrong=wrong)
        status = "completed"

    return {
        "target_ancestry": target,
        "wrong_ancestry": wrong,
        "runner_wrong_label": runner_wrong,
        "status": status,
        "out_dir": str(pair_out_dir),
        "aliased_manifest": str(aliased_manifest),
        "command": " ".join(command),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-dir", type=Path, default=DEFAULT_PANEL_DIR)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument(
        "--shared-reference-panel",
        action="store_true",
        help="Use one R=100 reference panel for LD/Saleem/pairwise and template rows (disjoint from beacon only).",
    )
    parser.add_argument("--target-ancestries", default="EUR,AFR,EAS,SAS,AMR")
    parser.add_argument("--wrong-ancestries", default="EUR,AFR,EAS,SAS,AMR")
    parser.add_argument("--beacon-sizes", default="25,50,100,250")
    parser.add_argument("--max-windows", type=int, default=6)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--matched-fractions", default="0,25,50,75,100")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--proposals-per-epoch", type=int, default=5000)
    parser.add_argument("--early-stopping", type=int, default=5)
    parser.add_argument("--dry-run", action="store_true", help="Print and record planned comparisons without launching runs.")
    parser.add_argument("--overwrite", action="store_true", help="Re-run comparisons whose output metrics already exist.")
    add_saleem_arguments(parser)
    args = parser.parse_args()

    args.out_dir = args.out_dir or (DEFAULT_SHARED_OUT_DIR if args.shared_reference_panel else DEFAULT_OUT_DIR)

    source_rows = _read_csv(args.panel_dir / "ancestry_panel_manifest.csv")
    available = {row["ancestry"].upper() for row in source_rows}
    target_ancestries = _split_csv(args.target_ancestries)
    wrong_ancestries = _split_csv(args.wrong_ancestries)
    missing = (set(target_ancestries) | set(wrong_ancestries)) - available
    if missing:
        raise ValueError(f"Requested ancestries are missing from {args.panel_dir}: {sorted(missing)}")

    pairs = _comparison_pairs(target_ancestries, wrong_ancestries)
    summary_rows = []
    for target, wrong in pairs:
        print(f"launching {target} target with {wrong} wrong reference")
        summary_rows.append(_run_comparison(args=args, source_rows=source_rows, target=target, wrong=wrong))

    if not args.dry_run:
        combined_path = _rebuild_combined_metrics(args=args, pairs=pairs)
        if combined_path is not None:
            print(f"wrote {combined_path}")

    if summary_rows:
        suffix = "dry_run" if args.dry_run else "summary"
        summary_path = args.out_dir / f"multiancestry_all_wrong_references_{suffix}.csv"
        _write_csv(summary_path, summary_rows)
        print(f"wrote {summary_path}")
    else:
        print("No target/wrong comparisons requested.")


if __name__ == "__main__":
    main()
