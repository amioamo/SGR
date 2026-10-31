# Partial-Genotype Linkage on Reconstructed Rows

## Configuration

- windows: 6
- replicates: 3
- beacon sizes: 25,50,100,250
- known SNP counts: 10,25,50,100
- methods: matched_template_n100

## Aggregate Linkage Accuracy

| N | method | known SNPs | runs | top-1 | top-5 | mean rank |
|---:|---|---:|---:|---:|---:|---:|
| 25 | independent_count | 10 | 18 | 0.0444 | 0.2444 | 11.56 |
| 25 | independent_count | 25 | 18 | 0.0622 | 0.2533 | 11.73 |
| 25 | independent_count | 50 | 18 | 0.0844 | 0.3311 | 10.64 |
| 25 | independent_count | 100 | 18 | 0.1067 | 0.3711 | 9.38 |
| 25 | matched_template_n100 | 10 | 18 | 0.1667 | 0.5178 | 7.24 |
| 25 | matched_template_n100 | 25 | 18 | 0.2133 | 0.5867 | 5.95 |
| 25 | matched_template_n100 | 50 | 18 | 0.2578 | 0.7044 | 4.88 |
| 25 | matched_template_n100 | 100 | 18 | 0.2956 | 0.7578 | 4.14 |
| 25 | pairwise_refinement | 10 | 18 | 0.1733 | 0.4756 | 7.64 |
| 25 | pairwise_refinement | 25 | 18 | 0.2156 | 0.5333 | 6.66 |
| 25 | pairwise_refinement | 50 | 18 | 0.2378 | 0.6644 | 5.40 |
| 25 | pairwise_refinement | 100 | 18 | 0.3000 | 0.7333 | 4.45 |
| 25 | saleem_reconstruction_attack | 10 | 18 | 0.1222 | 0.4422 | 8.26 |
| 25 | saleem_reconstruction_attack | 25 | 18 | 0.2000 | 0.5667 | 6.31 |
| 25 | saleem_reconstruction_attack | 50 | 18 | 0.2133 | 0.6089 | 5.87 |
| 25 | saleem_reconstruction_attack | 100 | 18 | 0.2711 | 0.6600 | 5.19 |
| 25 | shuffled_template_n100 | 10 | 18 | 0.0489 | 0.2333 | 12.14 |
| 25 | shuffled_template_n100 | 25 | 18 | 0.0822 | 0.2756 | 11.62 |
| 25 | shuffled_template_n100 | 50 | 18 | 0.0778 | 0.2956 | 10.70 |
| 25 | shuffled_template_n100 | 100 | 18 | 0.0889 | 0.3867 | 9.10 |
| 50 | independent_count | 10 | 18 | 0.0278 | 0.1322 | 22.78 |
| 50 | independent_count | 25 | 18 | 0.0489 | 0.1467 | 22.44 |
| 50 | independent_count | 50 | 18 | 0.0411 | 0.1778 | 20.07 |
| 50 | independent_count | 100 | 18 | 0.0656 | 0.2489 | 16.69 |
| 50 | matched_template_n100 | 10 | 18 | 0.0956 | 0.3400 | 13.68 |
| 50 | matched_template_n100 | 25 | 18 | 0.1478 | 0.4578 | 10.28 |
| 50 | matched_template_n100 | 50 | 18 | 0.1911 | 0.5356 | 8.48 |
| 50 | matched_template_n100 | 100 | 18 | 0.2500 | 0.6456 | 6.38 |
| 50 | pairwise_refinement | 10 | 18 | 0.0811 | 0.2878 | 14.96 |
| 50 | pairwise_refinement | 25 | 18 | 0.1111 | 0.4000 | 12.19 |
| 50 | pairwise_refinement | 50 | 18 | 0.1578 | 0.4800 | 9.72 |
| 50 | pairwise_refinement | 100 | 18 | 0.2200 | 0.5789 | 7.84 |
| 50 | saleem_reconstruction_attack | 10 | 18 | 0.0811 | 0.3000 | 14.67 |
| 50 | saleem_reconstruction_attack | 25 | 18 | 0.1278 | 0.3911 | 11.93 |
| 50 | saleem_reconstruction_attack | 50 | 18 | 0.1756 | 0.4778 | 9.87 |
| 50 | saleem_reconstruction_attack | 100 | 18 | 0.2022 | 0.5444 | 8.52 |
| 50 | shuffled_template_n100 | 10 | 18 | 0.0333 | 0.1222 | 23.13 |
| 50 | shuffled_template_n100 | 25 | 18 | 0.0378 | 0.1656 | 21.65 |
| 50 | shuffled_template_n100 | 50 | 18 | 0.0467 | 0.2011 | 19.65 |
| 50 | shuffled_template_n100 | 100 | 18 | 0.0733 | 0.2211 | 18.32 |
| 100 | independent_count | 10 | 18 | 0.0156 | 0.0844 | 45.32 |
| 100 | independent_count | 25 | 18 | 0.0183 | 0.0817 | 41.83 |
| 100 | independent_count | 50 | 18 | 0.0333 | 0.1117 | 37.94 |
| 100 | independent_count | 100 | 18 | 0.0367 | 0.1472 | 34.05 |
| 100 | matched_template_n100 | 10 | 18 | 0.0550 | 0.2233 | 25.96 |
| 100 | matched_template_n100 | 25 | 18 | 0.0989 | 0.3250 | 19.07 |
| 100 | matched_template_n100 | 50 | 18 | 0.1317 | 0.4122 | 14.54 |
| 100 | matched_template_n100 | 100 | 18 | 0.1678 | 0.5061 | 11.35 |
| 100 | pairwise_refinement | 10 | 18 | 0.0467 | 0.1656 | 31.05 |
| 100 | pairwise_refinement | 25 | 18 | 0.0633 | 0.2417 | 24.62 |
| 100 | pairwise_refinement | 50 | 18 | 0.0944 | 0.3200 | 19.66 |
| 100 | pairwise_refinement | 100 | 18 | 0.1300 | 0.3911 | 16.35 |
| 100 | saleem_reconstruction_attack | 10 | 18 | 0.0556 | 0.1900 | 29.59 |
| 100 | saleem_reconstruction_attack | 25 | 18 | 0.0867 | 0.2706 | 23.01 |
| 100 | saleem_reconstruction_attack | 50 | 18 | 0.1122 | 0.3456 | 18.90 |
| 100 | saleem_reconstruction_attack | 100 | 18 | 0.1294 | 0.3889 | 16.00 |
| 100 | shuffled_template_n100 | 10 | 18 | 0.0178 | 0.0667 | 45.34 |
| 100 | shuffled_template_n100 | 25 | 18 | 0.0222 | 0.0978 | 42.11 |
| 100 | shuffled_template_n100 | 50 | 18 | 0.0250 | 0.1017 | 38.85 |
| 100 | shuffled_template_n100 | 100 | 18 | 0.0372 | 0.1417 | 34.66 |
| 250 | independent_count | 10 | 18 | 0.0087 | 0.0360 | 110.74 |
| 250 | independent_count | 25 | 18 | 0.0122 | 0.0462 | 103.96 |
| 250 | independent_count | 50 | 18 | 0.0178 | 0.0589 | 94.25 |
| 250 | independent_count | 100 | 18 | 0.0233 | 0.0858 | 82.27 |
| 250 | matched_template_n100 | 10 | 18 | 0.0218 | 0.1078 | 62.93 |
| 250 | matched_template_n100 | 25 | 18 | 0.0549 | 0.1871 | 44.95 |
| 250 | matched_template_n100 | 50 | 18 | 0.0718 | 0.2658 | 33.10 |
| 250 | matched_template_n100 | 100 | 18 | 0.1107 | 0.3513 | 24.72 |
| 250 | pairwise_refinement | 10 | 18 | 0.0162 | 0.0762 | 83.68 |
| 250 | pairwise_refinement | 25 | 18 | 0.0293 | 0.1158 | 67.56 |
| 250 | pairwise_refinement | 50 | 18 | 0.0456 | 0.1518 | 55.48 |
| 250 | pairwise_refinement | 100 | 18 | 0.0538 | 0.2000 | 46.03 |
| 250 | saleem_reconstruction_attack | 10 | 18 | 0.0247 | 0.0898 | 72.50 |
| 250 | saleem_reconstruction_attack | 25 | 18 | 0.0418 | 0.1496 | 55.14 |
| 250 | saleem_reconstruction_attack | 50 | 18 | 0.0671 | 0.1971 | 45.80 |
| 250 | saleem_reconstruction_attack | 100 | 18 | 0.0802 | 0.2507 | 38.80 |
| 250 | shuffled_template_n100 | 10 | 18 | 0.0071 | 0.0378 | 109.87 |
| 250 | shuffled_template_n100 | 25 | 18 | 0.0109 | 0.0464 | 102.69 |
| 250 | shuffled_template_n100 | 50 | 18 | 0.0127 | 0.0507 | 95.19 |
| 250 | shuffled_template_n100 | 100 | 18 | 0.0218 | 0.0802 | 82.30 |

## Interpretation

Each participant is paired to a reconstructed row after anonymous Hungarian alignment. An attacker with carrier status at k randomly selected SNPs ranks reconstructed rows by partial Hamming distance and is scored on whether the aligned row is recovered in the top ranks.
