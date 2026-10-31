# Stage Ablation On The Matched Reference Panel

## Configuration

- windows: 6
- replicates: 3
- beacon sizes: 25,50,100,250
- panel size: 100
- reference panel mode: shared

Stage 1 of the full method is greedy template selection. Stage 2 is exact count projection.
`reference_only` samples N panel rows and does not use y.
`random_template_projection` uses that same random draw as the template, then applies Stage 2.
`greedy_template_without_projection` keeps the Stage 1 template and skips Stage 2.
Count violations are expected for the two methods that skip projection; their F1 is still reported.

## Aggregate Metrics

| N | method | n_runs | mean_f1 | median_f1 | mean_hamming | mean_count_mae | mean_max_count_error | invalid |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 25 | reference_only | 18 | 0.712430 | 0.694808 | 885.83 | 1.2372 | 6.1111 | 18 |
| 25 | random_template_projection | 18 | 0.674592 | 0.644406 | 1027.22 | 0.0000 | 0.0000 | 0 |
| 25 | greedy_template_without_projection | 18 | 0.796554 | 0.746017 | 577.33 | 0.3124 | 3.5000 | 18 |
| 25 | template_count_synthesis | 18 | 0.788542 | 0.731127 | 612.78 | 0.0000 | 0.0000 | 0 |
| 50 | reference_only | 18 | 0.721195 | 0.710898 | 1623.56 | 2.0793 | 8.2222 | 18 |
| 50 | random_template_projection | 18 | 0.669061 | 0.657173 | 1877.56 | 0.0000 | 0.0000 | 0 |
| 50 | greedy_template_without_projection | 18 | 0.816740 | 0.772217 | 1056.61 | 0.4323 | 4.6111 | 18 |
| 50 | template_count_synthesis | 18 | 0.809228 | 0.767371 | 1101.89 | 0.0000 | 0.0000 | 0 |
| 100 | reference_only | 18 | 0.770216 | 0.774035 | 2643.67 | 3.0651 | 12.1667 | 18 |
| 100 | random_template_projection | 18 | 0.717731 | 0.730072 | 3157.44 | 0.0000 | 0.0000 | 0 |
| 100 | greedy_template_without_projection | 18 | 0.832127 | 0.797198 | 1924.11 | 0.6622 | 7.6111 | 18 |
| 100 | template_count_synthesis | 18 | 0.824958 | 0.791479 | 2003.11 | 0.0000 | 0.0000 | 0 |
| 250 | reference_only | 18 | 0.771340 | 0.758113 | 6535.39 | 7.5119 | 29.9444 | 18 |
| 250 | random_template_projection | 18 | 0.716465 | 0.716122 | 8103.33 | 0.0000 | 0.0000 | 0 |
| 250 | greedy_template_without_projection | 18 | 0.841549 | 0.800224 | 4579.28 | 1.3592 | 14.1111 | 18 |
| 250 | template_count_synthesis | 18 | 0.835373 | 0.795644 | 4761.00 | 0.0000 | 0.0000 | 0 |

## Gain Versus Each Baseline

Positive `full_minus_baseline_f1` is the F1 gained by using both stages instead of that baseline.
`greedy_template_without_projection` isolates the Stage 2 projection gain.
`random_template_projection` isolates the Stage 1 greedy-selection gain.

| N | baseline | mean_baseline_f1 | mean_full_f1 | full_minus_baseline_f1 | wins_for_full | ties | losses |
|---:|---|---:|---:|---:|---:|---:|---:|
| 25 | reference_only | 0.712430 | 0.788542 | 0.076112 | 16 | 0 | 2 |
| 25 | random_template_projection | 0.674592 | 0.788542 | 0.113950 | 18 | 0 | 0 |
| 25 | greedy_template_without_projection | 0.796554 | 0.788542 | -0.008012 | 3 | 0 | 15 |
| 50 | reference_only | 0.721195 | 0.809228 | 0.088032 | 17 | 0 | 1 |
| 50 | random_template_projection | 0.669061 | 0.809228 | 0.140167 | 18 | 0 | 0 |
| 50 | greedy_template_without_projection | 0.816740 | 0.809228 | -0.007512 | 0 | 0 | 18 |
| 100 | reference_only | 0.770216 | 0.824958 | 0.054742 | 16 | 0 | 2 |
| 100 | random_template_projection | 0.717731 | 0.824958 | 0.107227 | 18 | 0 | 0 |
| 100 | greedy_template_without_projection | 0.832127 | 0.824958 | -0.007169 | 1 | 0 | 17 |
| 250 | reference_only | 0.771340 | 0.835373 | 0.064033 | 18 | 0 | 0 |
| 250 | random_template_projection | 0.716465 | 0.835373 | 0.118908 | 18 | 0 | 0 |
| 250 | greedy_template_without_projection | 0.841549 | 0.835373 | -0.006176 | 1 | 0 | 17 |

## Participant Similarity

Similarity is `1 - Hamming / m` after the same one-to-one Hungarian alignment used for carrier F1.
`sim>=0.95` is the fraction of aligned participants at or above 0.95.

| N | method | n_runs | participants | mean_similarity | sim>=0.95 |
|---:|---|---:|---:|---:|---:|
| 25 | reference_only | 18 | 450 | 0.929133 | 0.533333 |
| 25 | random_template_projection | 18 | 450 | 0.917822 | 0.377778 |
| 25 | greedy_template_without_projection | 18 | 450 | 0.953813 | 0.653333 |
| 25 | template_count_synthesis | 18 | 450 | 0.950978 | 0.635556 |
| 50 | reference_only | 18 | 900 | 0.935058 | 0.611111 |
| 50 | random_template_projection | 18 | 900 | 0.924898 | 0.431111 |
| 50 | greedy_template_without_projection | 18 | 900 | 0.957736 | 0.677778 |
| 50 | template_count_synthesis | 18 | 900 | 0.955924 | 0.670000 |
| 100 | reference_only | 18 | 1800 | 0.947127 | 0.684444 |
| 100 | random_template_projection | 18 | 1800 | 0.936851 | 0.538889 |
| 100 | greedy_template_without_projection | 18 | 1800 | 0.961518 | 0.721111 |
| 100 | template_count_synthesis | 18 | 1800 | 0.959938 | 0.706111 |
| 250 | reference_only | 18 | 4500 | 0.947717 | 0.689556 |
| 250 | random_template_projection | 18 | 4500 | 0.935173 | 0.525333 |
| 250 | greedy_template_without_projection | 18 | 4500 | 0.963366 | 0.737111 |
| 250 | template_count_synthesis | 18 | 4500 | 0.961912 | 0.726889 |

## Paired Inference

Each pair is one window by replicate. The 95% CI is a percentile bootstrap of the mean paired difference (10000 resamples, seed 0). The sign test is two-sided and drops ties. The Wilcoxon signed-rank test is two-sided and drops exact zeros.
For F1 and sim>=0.95, delta is full method minus baseline, and a win means the full method is higher.
Count-violation runs are retained.

| N | comparison | metric | n | mean_delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |
|---:|---|---|---:|---:|---|---:|---:|---:|
| 25 | template_count_synthesis_minus_greedy_template_without_projection | f1 | 18 | -0.008012 | [-0.012664, -0.003736] | 3/0/15 | 0.007538 | 0.001579 |
| 25 | template_count_synthesis_minus_random_template_projection | f1 | 18 | 0.113950 | [0.078801, 0.153241] | 18/0/0 | 7.629e-06 | 7.629e-06 |
| 25 | template_count_synthesis_minus_reference_only | f1 | 18 | 0.076112 | [0.043623, 0.112851] | 16/0/2 | 0.001312 | 7.629e-05 |
| 25 | template_count_synthesis_minus_greedy_template_without_projection | proportion_similarity_ge_095 | 18 | -0.017778 | [-0.042222, 0.004444] | 2/10/6 | 0.2891 | 0.1797 |
| 25 | template_count_synthesis_minus_random_template_projection | proportion_similarity_ge_095 | 18 | 0.257778 | [0.171111, 0.355556] | 16/2/0 | 3.052e-05 | 0.0004349 |
| 25 | template_count_synthesis_minus_reference_only | proportion_similarity_ge_095 | 18 | 0.102222 | [0.055556, 0.153333] | 13/4/1 | 0.001831 | 0.001198 |
| 50 | template_count_synthesis_minus_greedy_template_without_projection | f1 | 18 | -0.007512 | [-0.010644, -0.004837] | 0/0/18 | 7.629e-06 | 7.629e-06 |
| 50 | template_count_synthesis_minus_random_template_projection | f1 | 18 | 0.140167 | [0.093797, 0.196730] | 18/0/0 | 7.629e-06 | 7.629e-06 |
| 50 | template_count_synthesis_minus_reference_only | f1 | 18 | 0.088032 | [0.050625, 0.132029] | 17/0/1 | 0.000145 | 3.815e-05 |
| 50 | template_count_synthesis_minus_greedy_template_without_projection | proportion_similarity_ge_095 | 18 | -0.007778 | [-0.021111, 0.003333] | 4/10/4 | 1 | 0.3047 |
| 50 | template_count_synthesis_minus_random_template_projection | proportion_similarity_ge_095 | 18 | 0.238889 | [0.164444, 0.318889] | 17/1/0 | 1.526e-05 | 1.526e-05 |
| 50 | template_count_synthesis_minus_reference_only | proportion_similarity_ge_095 | 18 | 0.058889 | [0.030000, 0.087778] | 15/1/2 | 0.00235 | 0.002085 |
| 100 | template_count_synthesis_minus_greedy_template_without_projection | f1 | 18 | -0.007169 | [-0.010044, -0.004696] | 1/0/17 | 0.000145 | 2.289e-05 |
| 100 | template_count_synthesis_minus_random_template_projection | f1 | 18 | 0.107227 | [0.067070, 0.163462] | 18/0/0 | 7.629e-06 | 7.629e-06 |
| 100 | template_count_synthesis_minus_reference_only | f1 | 18 | 0.054742 | [0.020561, 0.104960] | 16/0/2 | 0.001312 | 0.0003281 |
| 100 | template_count_synthesis_minus_greedy_template_without_projection | proportion_similarity_ge_095 | 18 | -0.015000 | [-0.027778, -0.005000] | 2/8/8 | 0.1094 | 0.01562 |
| 100 | template_count_synthesis_minus_random_template_projection | proportion_similarity_ge_095 | 18 | 0.167222 | [0.120556, 0.215556] | 17/1/0 | 1.526e-05 | 1.526e-05 |
| 100 | template_count_synthesis_minus_reference_only | proportion_similarity_ge_095 | 18 | 0.021667 | [0.000000, 0.043889] | 11/1/6 | 0.3323 | 0.1408 |
| 250 | template_count_synthesis_minus_greedy_template_without_projection | f1 | 18 | -0.006176 | [-0.008317, -0.004147] | 1/0/17 | 0.000145 | 1.526e-05 |
| 250 | template_count_synthesis_minus_random_template_projection | f1 | 18 | 0.118908 | [0.088628, 0.154714] | 18/0/0 | 7.629e-06 | 7.629e-06 |
| 250 | template_count_synthesis_minus_reference_only | f1 | 18 | 0.064033 | [0.042963, 0.088621] | 18/0/0 | 7.629e-06 | 7.629e-06 |
| 250 | template_count_synthesis_minus_greedy_template_without_projection | proportion_similarity_ge_095 | 18 | -0.010222 | [-0.016889, -0.003778] | 3/3/12 | 0.03516 | 0.01231 |
| 250 | template_count_synthesis_minus_random_template_projection | proportion_similarity_ge_095 | 18 | 0.201556 | [0.138222, 0.274672] | 17/0/1 | 0.000145 | 0.0002328 |
| 250 | template_count_synthesis_minus_reference_only | proportion_similarity_ge_095 | 18 | 0.037333 | [0.018222, 0.056444] | 13/1/4 | 0.04904 | 0.006482 |

## Projection Edit Fraction

For a template `T`, the edit fraction is `||y - c(T)||_1 / (N m)`.
Stage 2 changes exactly that fraction of entries, because each column flips `|y_j - c_j(T)|` bits.
The random template is the `reference_only` sample. The greedy template is Stage 1.
Delta is greedy minus random, so a negative value means the greedy template needs fewer edits.
A win means the greedy template has the smaller edit fraction.

| N | template | mean_edit_fraction |
|---:|---|---:|
| 25 | reference_only | 0.049489 |
| 25 | greedy_template_without_projection | 0.012498 |
| 50 | reference_only | 0.041587 |
| 50 | greedy_template_without_projection | 0.008647 |
| 100 | reference_only | 0.030651 |
| 100 | greedy_template_without_projection | 0.006622 |
| 250 | reference_only | 0.030048 |
| 250 | greedy_template_without_projection | 0.005437 |

| N | comparison | n | mean_delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |
|---:|---|---:|---:|---|---:|---:|---:|
| 25 | greedy_template_without_projection_minus_reference_only | 18 | -0.036991 | [-0.048840, -0.026320] | 18/0/0 | 7.629e-06 | 7.629e-06 |
| 50 | greedy_template_without_projection_minus_reference_only | 18 | -0.032940 | [-0.040985, -0.025800] | 18/0/0 | 7.629e-06 | 7.629e-06 |
| 100 | greedy_template_without_projection_minus_reference_only | 18 | -0.024029 | [-0.030947, -0.018079] | 18/0/0 | 7.629e-06 | 7.629e-06 |
| 250 | greedy_template_without_projection_minus_reference_only | 18 | -0.024611 | [-0.029182, -0.020940] | 18/0/0 | 7.629e-06 | 7.629e-06 |

## Claim Boundary

Rows are aligned with Hungarian matching before F1 is computed, including runs whose column sums
do not match y. `invalid_run` records that count violation and does not drop the F1.
Reference panel mode: `shared`.

