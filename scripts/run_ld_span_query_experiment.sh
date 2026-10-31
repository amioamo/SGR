#!/usr/bin/env bash
# Additional regions x dense / LD strength / sparse / long-range.
# Fast methods write the metrics file; Saleem is merged in a second pass.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
python scripts/prepare_ld_span_query_panels.py \
  --chroms 1,20,22 \
  --regions-per-chrom 2 \
  --snp-count 200 \
  --out-dir data/ld_span_queries
OUT=results/block_inheritance_v1/ld_span_queries_shared_panel
python scripts/run_ld_span_query_evaluation.py \
  --manifest data/ld_span_queries/ld_span_query_manifest.csv \
  --max-replicates 3 \
  --beacon-sizes 25,50 \
  --no-run-saleem \
  --out-dir "$OUT"
python scripts/run_ld_span_query_evaluation.py \
  --manifest data/ld_span_queries/ld_span_query_manifest.csv \
  --max-replicates 3 \
  --beacon-sizes 25,50 \
  --saleem-only \
  --merge-metrics "$OUT/ld_span_query_metrics.csv" \
  --out-dir "$OUT"
