# LD, Density, And Genomic-Span Query Evaluation

Fixed SNP count across additional genomic regions. `dense` is the short-span, high-density reference. `ld_r2_080`, `ld_r2_050`, and `ld_r2_010` cap pairwise carrier r^2 at 0.80, 0.50, and 0.10. `sparse` enforces a 50 kb minimum gap. `long_range` spaces the same SNP count across 10 Mb.

## Configuration

- query sets: 36
- replicates: 3
- beacon sizes: [25, 50]
- reference panel: shared
- template selection: proportional_greedy_with_replacement

## Realized Panel Structure And Carrier-Level F1

| N | regime | method | n | span Mb | gap kb | SNPs/Mb | mean |r| | max r2 | mean F1 | invalid |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 25 | dense | independent_count | 18 | 0.085 | 0.242 | 2664.8 | 0.229 | 1.000 | 0.379 | 18 |
| 25 | dense | matched_template_n100 | 18 | 0.085 | 0.242 | 2664.8 | 0.229 | 1.000 | 0.792 | 18 |
| 25 | dense | pairwise_refinement | 18 | 0.085 | 0.242 | 2664.8 | 0.229 | 1.000 | 0.746 | 18 |
| 25 | dense | saleem_reconstruction_attack | 18 | 0.085 | 0.242 | 2664.8 | 0.229 | 1.000 | 0.705 | 0 |
| 25 | dense | shuffled_template_n100 | 18 | 0.085 | 0.242 | 2664.8 | 0.229 | 1.000 | 0.379 | 18 |
| 25 | ld_r2_010 | independent_count | 18 | 1.130 | 3.013 | 211.1 | 0.045 | 0.100 | 0.401 | 18 |
| 25 | ld_r2_010 | matched_template_n100 | 18 | 1.130 | 3.013 | 211.1 | 0.045 | 0.100 | 0.433 | 18 |
| 25 | ld_r2_010 | pairwise_refinement | 18 | 1.130 | 3.013 | 211.1 | 0.045 | 0.100 | 0.423 | 18 |
| 25 | ld_r2_010 | saleem_reconstruction_attack | 18 | 1.130 | 3.013 | 211.1 | 0.045 | 0.100 | 0.427 | 0 |
| 25 | ld_r2_010 | shuffled_template_n100 | 18 | 1.130 | 3.013 | 211.1 | 0.045 | 0.100 | 0.402 | 18 |
| 25 | ld_r2_050 | independent_count | 18 | 0.536 | 1.438 | 443.0 | 0.061 | 0.498 | 0.383 | 18 |
| 25 | ld_r2_050 | matched_template_n100 | 18 | 0.536 | 1.438 | 443.0 | 0.061 | 0.498 | 0.504 | 18 |
| 25 | ld_r2_050 | pairwise_refinement | 18 | 0.536 | 1.438 | 443.0 | 0.061 | 0.498 | 0.486 | 18 |
| 25 | ld_r2_050 | saleem_reconstruction_attack | 18 | 0.536 | 1.438 | 443.0 | 0.061 | 0.498 | 0.473 | 0 |
| 25 | ld_r2_050 | shuffled_template_n100 | 18 | 0.536 | 1.438 | 443.0 | 0.061 | 0.498 | 0.385 | 18 |
| 25 | ld_r2_080 | independent_count | 18 | 0.356 | 0.989 | 630.7 | 0.073 | 0.798 | 0.383 | 18 |
| 25 | ld_r2_080 | matched_template_n100 | 18 | 0.356 | 0.989 | 630.7 | 0.073 | 0.798 | 0.552 | 18 |
| 25 | ld_r2_080 | pairwise_refinement | 18 | 0.356 | 0.989 | 630.7 | 0.073 | 0.798 | 0.532 | 18 |
| 25 | ld_r2_080 | saleem_reconstruction_attack | 18 | 0.356 | 0.989 | 630.7 | 0.073 | 0.798 | 0.520 | 0 |
| 25 | ld_r2_080 | shuffled_template_n100 | 18 | 0.356 | 0.989 | 630.7 | 0.073 | 0.798 | 0.378 | 18 |
| 25 | long_range | independent_count | 18 | 9.999 | 46.418 | 20.0 | 0.042 | 1.000 | 0.399 | 18 |
| 25 | long_range | matched_template_n100 | 18 | 9.999 | 46.418 | 20.0 | 0.042 | 1.000 | 0.419 | 18 |
| 25 | long_range | pairwise_refinement | 18 | 9.999 | 46.418 | 20.0 | 0.042 | 1.000 | 0.415 | 18 |
| 25 | long_range | saleem_reconstruction_attack | 18 | 9.999 | 46.418 | 20.0 | 0.042 | 1.000 | 0.418 | 0 |
| 25 | long_range | shuffled_template_n100 | 18 | 9.999 | 46.418 | 20.0 | 0.042 | 1.000 | 0.400 | 18 |
| 25 | sparse | independent_count | 18 | 10.064 | 50.323 | 19.9 | 0.042 | 0.992 | 0.401 | 18 |
| 25 | sparse | matched_template_n100 | 18 | 10.064 | 50.323 | 19.9 | 0.042 | 0.992 | 0.420 | 18 |
| 25 | sparse | pairwise_refinement | 18 | 10.064 | 50.323 | 19.9 | 0.042 | 0.992 | 0.418 | 18 |
| 25 | sparse | saleem_reconstruction_attack | 18 | 10.064 | 50.323 | 19.9 | 0.042 | 0.992 | 0.419 | 0 |
| 25 | sparse | shuffled_template_n100 | 18 | 10.064 | 50.323 | 19.9 | 0.042 | 0.992 | 0.400 | 18 |
| 50 | dense | independent_count | 18 | 0.085 | 0.242 | 2664.8 | 0.229 | 1.000 | 0.375 | 18 |
| 50 | dense | matched_template_n100 | 18 | 0.085 | 0.242 | 2664.8 | 0.229 | 1.000 | 0.828 | 18 |
| 50 | dense | pairwise_refinement | 18 | 0.085 | 0.242 | 2664.8 | 0.229 | 1.000 | 0.764 | 18 |
| 50 | dense | saleem_reconstruction_attack | 18 | 0.085 | 0.242 | 2664.8 | 0.229 | 1.000 | 0.710 | 0 |
| 50 | dense | shuffled_template_n100 | 18 | 0.085 | 0.242 | 2664.8 | 0.229 | 1.000 | 0.374 | 18 |
| 50 | ld_r2_010 | independent_count | 18 | 1.130 | 3.013 | 211.1 | 0.045 | 0.100 | 0.405 | 18 |
| 50 | ld_r2_010 | matched_template_n100 | 18 | 1.130 | 3.013 | 211.1 | 0.045 | 0.100 | 0.437 | 18 |
| 50 | ld_r2_010 | pairwise_refinement | 18 | 1.130 | 3.013 | 211.1 | 0.045 | 0.100 | 0.425 | 18 |
| 50 | ld_r2_010 | saleem_reconstruction_attack | 18 | 1.130 | 3.013 | 211.1 | 0.045 | 0.100 | 0.423 | 0 |
| 50 | ld_r2_010 | shuffled_template_n100 | 18 | 1.130 | 3.013 | 211.1 | 0.045 | 0.100 | 0.400 | 18 |
| 50 | ld_r2_050 | independent_count | 18 | 0.536 | 1.438 | 443.0 | 0.061 | 0.498 | 0.386 | 18 |
| 50 | ld_r2_050 | matched_template_n100 | 18 | 0.536 | 1.438 | 443.0 | 0.061 | 0.498 | 0.515 | 18 |
| 50 | ld_r2_050 | pairwise_refinement | 18 | 0.536 | 1.438 | 443.0 | 0.061 | 0.498 | 0.490 | 18 |
| 50 | ld_r2_050 | saleem_reconstruction_attack | 18 | 0.536 | 1.438 | 443.0 | 0.061 | 0.498 | 0.470 | 0 |
| 50 | ld_r2_050 | shuffled_template_n100 | 18 | 0.536 | 1.438 | 443.0 | 0.061 | 0.498 | 0.384 | 18 |
| 50 | ld_r2_080 | independent_count | 18 | 0.356 | 0.989 | 630.7 | 0.073 | 0.798 | 0.381 | 18 |
| 50 | ld_r2_080 | matched_template_n100 | 18 | 0.356 | 0.989 | 630.7 | 0.073 | 0.798 | 0.571 | 18 |
| 50 | ld_r2_080 | pairwise_refinement | 18 | 0.356 | 0.989 | 630.7 | 0.073 | 0.798 | 0.547 | 18 |
| 50 | ld_r2_080 | saleem_reconstruction_attack | 18 | 0.356 | 0.989 | 630.7 | 0.073 | 0.798 | 0.517 | 0 |
| 50 | ld_r2_080 | shuffled_template_n100 | 18 | 0.356 | 0.989 | 630.7 | 0.073 | 0.798 | 0.377 | 18 |
| 50 | long_range | independent_count | 18 | 9.999 | 46.418 | 20.0 | 0.042 | 1.000 | 0.401 | 18 |
| 50 | long_range | matched_template_n100 | 18 | 9.999 | 46.418 | 20.0 | 0.042 | 1.000 | 0.418 | 18 |
| 50 | long_range | pairwise_refinement | 18 | 9.999 | 46.418 | 20.0 | 0.042 | 1.000 | 0.417 | 18 |
| 50 | long_range | saleem_reconstruction_attack | 18 | 9.999 | 46.418 | 20.0 | 0.042 | 1.000 | 0.417 | 0 |
| 50 | long_range | shuffled_template_n100 | 18 | 9.999 | 46.418 | 20.0 | 0.042 | 1.000 | 0.399 | 18 |
| 50 | sparse | independent_count | 18 | 10.064 | 50.323 | 19.9 | 0.042 | 0.992 | 0.402 | 18 |
| 50 | sparse | matched_template_n100 | 18 | 10.064 | 50.323 | 19.9 | 0.042 | 0.992 | 0.420 | 18 |
| 50 | sparse | pairwise_refinement | 18 | 10.064 | 50.323 | 19.9 | 0.042 | 0.992 | 0.420 | 18 |
| 50 | sparse | saleem_reconstruction_attack | 18 | 10.064 | 50.323 | 19.9 | 0.042 | 0.992 | 0.421 | 0 |
| 50 | sparse | shuffled_template_n100 | 18 | 10.064 | 50.323 | 19.9 | 0.042 | 0.992 | 0.399 | 18 |

## Delta F1 Versus Local LD Refinement

| N | regime | method | mean delta | median delta | wins | ties | losses |
|---:|---|---|---:|---:|---:|---:|---:|
| 25 | dense | independent_count | -0.367 | -0.362 | 0 | 0 | 18 |
| 25 | dense | matched_template_n100 | 0.046 | 0.035 | 13 | 0 | 5 |
| 25 | dense | saleem_reconstruction_attack | -0.041 | -0.022 | 7 | 1 | 10 |
| 25 | dense | shuffled_template_n100 | -0.367 | -0.371 | 0 | 0 | 18 |
| 25 | ld_r2_010 | independent_count | -0.022 | -0.024 | 0 | 1 | 17 |
| 25 | ld_r2_010 | matched_template_n100 | 0.010 | 0.010 | 13 | 2 | 3 |
| 25 | ld_r2_010 | saleem_reconstruction_attack | 0.004 | 0.001 | 9 | 2 | 7 |
| 25 | ld_r2_010 | shuffled_template_n100 | -0.021 | -0.022 | 1 | 0 | 17 |
| 25 | ld_r2_050 | independent_count | -0.103 | -0.108 | 0 | 0 | 18 |
| 25 | ld_r2_050 | matched_template_n100 | 0.018 | 0.018 | 17 | 0 | 1 |
| 25 | ld_r2_050 | saleem_reconstruction_attack | -0.013 | -0.010 | 5 | 0 | 13 |
| 25 | ld_r2_050 | shuffled_template_n100 | -0.101 | -0.105 | 0 | 0 | 18 |
| 25 | ld_r2_080 | independent_count | -0.149 | -0.142 | 0 | 0 | 18 |
| 25 | ld_r2_080 | matched_template_n100 | 0.020 | 0.023 | 13 | 1 | 4 |
| 25 | ld_r2_080 | saleem_reconstruction_attack | -0.013 | -0.010 | 8 | 0 | 10 |
| 25 | ld_r2_080 | shuffled_template_n100 | -0.154 | -0.148 | 0 | 0 | 18 |
| 25 | long_range | independent_count | -0.016 | -0.016 | 1 | 0 | 17 |
| 25 | long_range | matched_template_n100 | 0.004 | 0.006 | 14 | 0 | 4 |
| 25 | long_range | saleem_reconstruction_attack | 0.003 | 0.000 | 8 | 3 | 7 |
| 25 | long_range | shuffled_template_n100 | -0.014 | -0.013 | 1 | 0 | 17 |
| 25 | sparse | independent_count | -0.017 | -0.017 | 1 | 1 | 16 |
| 25 | sparse | matched_template_n100 | 0.002 | 0.001 | 10 | 3 | 5 |
| 25 | sparse | saleem_reconstruction_attack | 0.001 | -0.001 | 8 | 0 | 10 |
| 25 | sparse | shuffled_template_n100 | -0.018 | -0.015 | 0 | 1 | 17 |
| 50 | dense | independent_count | -0.388 | -0.381 | 0 | 0 | 18 |
| 50 | dense | matched_template_n100 | 0.064 | 0.051 | 17 | 0 | 1 |
| 50 | dense | saleem_reconstruction_attack | -0.054 | -0.024 | 4 | 0 | 14 |
| 50 | dense | shuffled_template_n100 | -0.390 | -0.388 | 0 | 0 | 18 |
| 50 | ld_r2_010 | independent_count | -0.021 | -0.021 | 0 | 0 | 18 |
| 50 | ld_r2_010 | matched_template_n100 | 0.012 | 0.012 | 16 | 0 | 2 |
| 50 | ld_r2_010 | saleem_reconstruction_attack | -0.002 | -0.002 | 6 | 1 | 11 |
| 50 | ld_r2_010 | shuffled_template_n100 | -0.025 | -0.025 | 0 | 0 | 18 |
| 50 | ld_r2_050 | independent_count | -0.105 | -0.107 | 0 | 0 | 18 |
| 50 | ld_r2_050 | matched_template_n100 | 0.025 | 0.021 | 17 | 0 | 1 |
| 50 | ld_r2_050 | saleem_reconstruction_attack | -0.020 | -0.024 | 3 | 0 | 15 |
| 50 | ld_r2_050 | shuffled_template_n100 | -0.106 | -0.110 | 0 | 0 | 18 |
| 50 | ld_r2_080 | independent_count | -0.166 | -0.163 | 0 | 0 | 18 |
| 50 | ld_r2_080 | matched_template_n100 | 0.024 | 0.025 | 17 | 0 | 1 |
| 50 | ld_r2_080 | saleem_reconstruction_attack | -0.030 | -0.027 | 1 | 1 | 16 |
| 50 | ld_r2_080 | shuffled_template_n100 | -0.170 | -0.165 | 0 | 0 | 18 |
| 50 | long_range | independent_count | -0.016 | -0.017 | 0 | 0 | 18 |
| 50 | long_range | matched_template_n100 | 0.001 | 0.000 | 9 | 1 | 8 |
| 50 | long_range | saleem_reconstruction_attack | 0.000 | 0.000 | 9 | 2 | 7 |
| 50 | long_range | shuffled_template_n100 | -0.017 | -0.017 | 0 | 0 | 18 |
| 50 | sparse | independent_count | -0.018 | -0.017 | 2 | 0 | 16 |
| 50 | sparse | matched_template_n100 | 0.000 | 0.001 | 11 | 0 | 7 |
| 50 | sparse | saleem_reconstruction_attack | 0.001 | 0.002 | 11 | 0 | 7 |
| 50 | sparse | shuffled_template_n100 | -0.021 | -0.018 | 0 | 0 | 18 |

## Delta F1 Versus The Dense Panel In The Same Additional Region

| N | regime | method | mean delta | median delta | wins | ties | losses |
|---:|---|---|---:|---:|---:|---:|---:|
| 25 | ld_r2_010 | independent_count | 0.022 | 0.025 | 14 | 0 | 4 |
| 25 | ld_r2_010 | matched_template_n100 | -0.359 | -0.342 | 0 | 0 | 18 |
| 25 | ld_r2_010 | pairwise_refinement | -0.323 | -0.298 | 0 | 0 | 18 |
| 25 | ld_r2_010 | saleem_reconstruction_attack | -0.279 | -0.269 | 0 | 0 | 18 |
| 25 | ld_r2_010 | shuffled_template_n100 | 0.023 | 0.035 | 14 | 0 | 4 |
| 25 | ld_r2_050 | independent_count | 0.004 | 0.012 | 12 | 0 | 6 |
| 25 | ld_r2_050 | matched_template_n100 | -0.288 | -0.281 | 0 | 0 | 18 |
| 25 | ld_r2_050 | pairwise_refinement | -0.260 | -0.235 | 0 | 0 | 18 |
| 25 | ld_r2_050 | saleem_reconstruction_attack | -0.232 | -0.223 | 0 | 0 | 18 |
| 25 | ld_r2_050 | shuffled_template_n100 | 0.006 | 0.016 | 11 | 0 | 7 |
| 25 | ld_r2_080 | independent_count | 0.004 | 0.006 | 11 | 0 | 7 |
| 25 | ld_r2_080 | matched_template_n100 | -0.240 | -0.259 | 0 | 0 | 18 |
| 25 | ld_r2_080 | pairwise_refinement | -0.214 | -0.202 | 0 | 0 | 18 |
| 25 | ld_r2_080 | saleem_reconstruction_attack | -0.186 | -0.183 | 0 | 0 | 18 |
| 25 | ld_r2_080 | shuffled_template_n100 | -0.001 | 0.007 | 11 | 0 | 7 |
| 25 | long_range | independent_count | 0.020 | 0.016 | 12 | 0 | 6 |
| 25 | long_range | matched_template_n100 | -0.373 | -0.358 | 0 | 0 | 18 |
| 25 | long_range | pairwise_refinement | -0.331 | -0.329 | 0 | 0 | 18 |
| 25 | long_range | saleem_reconstruction_attack | -0.288 | -0.284 | 0 | 0 | 18 |
| 25 | long_range | shuffled_template_n100 | 0.022 | 0.024 | 12 | 0 | 6 |
| 25 | sparse | independent_count | 0.022 | 0.026 | 13 | 0 | 5 |
| 25 | sparse | matched_template_n100 | -0.372 | -0.366 | 0 | 0 | 18 |
| 25 | sparse | pairwise_refinement | -0.328 | -0.330 | 0 | 0 | 18 |
| 25 | sparse | saleem_reconstruction_attack | -0.287 | -0.283 | 0 | 0 | 18 |
| 25 | sparse | shuffled_template_n100 | 0.021 | 0.035 | 12 | 0 | 6 |
| 50 | ld_r2_010 | independent_count | 0.029 | 0.028 | 15 | 0 | 3 |
| 50 | ld_r2_010 | matched_template_n100 | -0.391 | -0.400 | 0 | 0 | 18 |
| 50 | ld_r2_010 | pairwise_refinement | -0.338 | -0.333 | 0 | 0 | 18 |
| 50 | ld_r2_010 | saleem_reconstruction_attack | -0.287 | -0.283 | 0 | 0 | 18 |
| 50 | ld_r2_010 | shuffled_template_n100 | 0.026 | 0.019 | 13 | 0 | 5 |
| 50 | ld_r2_050 | independent_count | 0.010 | 0.009 | 10 | 0 | 8 |
| 50 | ld_r2_050 | matched_template_n100 | -0.313 | -0.338 | 0 | 0 | 18 |
| 50 | ld_r2_050 | pairwise_refinement | -0.273 | -0.258 | 0 | 0 | 18 |
| 50 | ld_r2_050 | saleem_reconstruction_attack | -0.240 | -0.227 | 0 | 0 | 18 |
| 50 | ld_r2_050 | shuffled_template_n100 | 0.011 | 0.005 | 11 | 0 | 7 |
| 50 | ld_r2_080 | independent_count | 0.006 | 0.009 | 10 | 0 | 8 |
| 50 | ld_r2_080 | matched_template_n100 | -0.257 | -0.300 | 0 | 0 | 18 |
| 50 | ld_r2_080 | pairwise_refinement | -0.217 | -0.216 | 0 | 0 | 18 |
| 50 | ld_r2_080 | saleem_reconstruction_attack | -0.193 | -0.177 | 0 | 0 | 18 |
| 50 | ld_r2_080 | shuffled_template_n100 | 0.003 | 0.005 | 9 | 0 | 9 |
| 50 | long_range | independent_count | 0.025 | 0.026 | 13 | 0 | 5 |
| 50 | long_range | matched_template_n100 | -0.410 | -0.413 | 0 | 0 | 18 |
| 50 | long_range | pairwise_refinement | -0.347 | -0.341 | 0 | 0 | 18 |
| 50 | long_range | saleem_reconstruction_attack | -0.293 | -0.287 | 0 | 0 | 18 |
| 50 | long_range | shuffled_template_n100 | 0.025 | 0.027 | 12 | 0 | 6 |
| 50 | sparse | independent_count | 0.027 | 0.024 | 14 | 0 | 4 |
| 50 | sparse | matched_template_n100 | -0.408 | -0.413 | 0 | 0 | 18 |
| 50 | sparse | pairwise_refinement | -0.343 | -0.331 | 0 | 0 | 18 |
| 50 | sparse | saleem_reconstruction_attack | -0.289 | -0.273 | 0 | 0 | 18 |
| 50 | sparse | shuffled_template_n100 | 0.025 | 0.026 | 14 | 0 | 4 |

## Claim Boundary

These panels are additional regions on chromosomes 1, 20, and 22, disjoint from the six main dense windows. They vary query structure under a matched shared reference panel. They are not genome-wide random queries, empirical Beacon logs, or a replacement for the main dense-window analysis.

