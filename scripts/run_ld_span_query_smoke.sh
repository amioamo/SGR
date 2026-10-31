#!/usr/bin/env bash
# One additional region on chr22, short SNP panels, fast methods only.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
python scripts/prepare_ld_span_query_panels.py \
  --chroms 22 \
  --regions-per-chrom 1 \
  --snp-count 40 \
  --ld-pool-span-bp 8000000 \
  --out-dir data/ld_span_queries_smoke
python scripts/run_ld_span_query_evaluation.py \
  --manifest data/ld_span_queries_smoke/ld_span_query_manifest.csv \
  --max-replicates 1 \
  --beacon-sizes 25 \
  --no-run-saleem \
  --out-dir results/block_inheritance_v1/ld_span_queries_smoke
