# Exact Template Optimality At Small R And N

## Configuration

- windows: 6
- replicates: 3
- draws per window and replicate: 5
- settings R:N: 5:5,6:4,8:5,8:8,10:6,10:8
- reference panel: shared matched panel, disjoint from the beacon
- beacon rows: random subset of the N=25 cohort
- feasible templates: multisets of N rows from R reference rows, with replacement

`T*` minimizes `||y - c(T)||_1` over that family.
Absolute gap is `||y - c(T_greedy)||_1 - ||y - c(T*)||_1`.
Relative gap divides the absolute gap by `max(1, ||y - c(T*)||_1)`.
Optimal-hit rate is the fraction of instances with absolute gap 0.
The hit-rate interval is a percentile bootstrap of that fraction (10000 resamples, seed 0).

## Aggregate

| R | N | space | n | mean greedy L1 | mean optimal L1 | mean abs gap | median abs gap | max abs gap | mean rel gap | hit rate | hit 95% CI |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 5 | 5 | 126 | 90 | 171.43 | 170.09 | 1.34 | 0.00 | 32 | 0.0057 | 0.9000 | [0.8333, 0.9556] |
| 6 | 4 | 126 | 90 | 143.04 | 142.01 | 1.03 | 0.00 | 14 | 0.0157 | 0.8222 | [0.7444, 0.9000] |
| 8 | 5 | 792 | 90 | 145.64 | 143.30 | 2.34 | 0.00 | 41 | 0.0622 | 0.6667 | [0.5667, 0.7667] |
| 8 | 8 | 6435 | 90 | 235.06 | 233.29 | 1.77 | 0.00 | 27 | 0.0129 | 0.7333 | [0.6333, 0.8222] |
| 10 | 6 | 5005 | 90 | 125.63 | 122.68 | 2.96 | 0.00 | 32 | 0.0402 | 0.5444 | [0.4444, 0.6444] |
| 10 | 8 | 24310 | 90 | 205.18 | 202.76 | 2.42 | 0.00 | 22 | 0.0370 | 0.6222 | [0.5222, 0.7222] |

## Claim Boundary

This comparison is restricted to enumerated small panels. It evaluates the final count objective
of Stage 1 template selection and does not include the Stage 2 projection.

