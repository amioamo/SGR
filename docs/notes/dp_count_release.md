# Discrete-Laplace Count Release

Carrier counts are released with the discrete Laplace mechanism at sensitivity 1.
Noise is independent across SNPs and the released integer is clipped to [0, N].
Clipping is post-processing. The joint privacy claim uses basic composition,
not advanced composition and not a correlation-aware accountant.

Per-query accounting spends the nominal epsilon on every SNP, so m SNPs compose to m * epsilon.
Composed accounting splits the nominal epsilon evenly, so each SNP uses epsilon / m and the joint budget is epsilon.
Reconstruction is scored against the true carrier matrix. The attack only observes the clipped noisy counts.

## Configuration

- windows: 6
- replicates: 3
- beacon sizes: 25, 50, 100, 250
- nominal epsilons: 0.01, 0.1, 1
- methods: independent_count, pairwise_refinement, template_count_synthesis, saleem_reconstruction_attack
- reference panel: shared
- template selection: proportional_greedy_with_replacement
- pairwise epochs: 20
- pairwise proposals per epoch: 5000
- Saleem correlation/frequency epochs: 1000/500
- Saleem max outer iterations: 3000

## Aggregate

| N | profile | accounting | nominal ε | ε per count | composed ε | method | runs | F1 | SD | ΔF1 | count MAE | changed |
|---:|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 25 | exact_count |  |  |  |  | independent_count | 6 | 0.2474 | 0.1241 |  | 0.0000 | 0.0000 |
| 25 | exact_count |  |  |  |  | pairwise_refinement | 6 | 0.8981 | 0.0461 |  | 0.0000 | 0.0000 |
| 25 | exact_count |  |  |  |  | saleem_reconstruction_attack | 6 | 0.7303 | 0.1952 |  | 0.0000 | 0.0000 |
| 25 | exact_count |  |  |  |  | template_count_synthesis | 6 | 0.9515 | 0.0279 |  | 0.0000 | 0.0000 |
| 25 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | independent_count | 12 | 0.1981 | 0.0899 | -0.0201 | 11.8843 | 0.9172 |
| 25 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | pairwise_refinement | 12 | 0.2044 | 0.0945 | -0.6611 | 11.8843 | 0.9172 |
| 25 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | saleem_reconstruction_attack | 12 | 0.2030 | 0.0905 | -0.4963 | 11.8843 | 0.9172 |
| 25 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | template_count_synthesis | 12 | 0.2037 | 0.0935 | -0.7144 | 11.8843 | 0.9172 |
| 25 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | independent_count | 12 | 0.2281 | 0.0891 | -0.0123 | 5.5118 | 0.8790 |
| 25 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | pairwise_refinement | 12 | 0.3288 | 0.1249 | -0.5328 | 5.5118 | 0.8790 |
| 25 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | saleem_reconstruction_attack | 12 | 0.3138 | 0.0936 | -0.3952 | 5.5118 | 0.8790 |
| 25 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | template_count_synthesis | 12 | 0.3337 | 0.1331 | -0.5731 | 5.5118 | 0.8790 |
| 25 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | independent_count | 12 | 0.2723 | 0.0932 | 0.0000 | 0.7208 | 0.5030 |
| 25 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | pairwise_refinement | 12 | 0.6928 | 0.0975 | -0.1285 | 0.7208 | 0.5030 |
| 25 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | saleem_reconstruction_attack | 11 | 0.6216 | 0.1099 | -0.1089 | 0.7240 | 0.5069 |
| 25 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | template_count_synthesis | 12 | 0.7341 | 0.1246 | -0.1324 | 0.7208 | 0.5030 |
| 25 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | independent_count | 6 | 0.2474 | 0.1241 | 0.0000 | 0.0000 | 0.0000 |
| 25 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | pairwise_refinement | 6 | 0.8981 | 0.0461 | 0.0000 | 0.0000 | 0.0000 |
| 25 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | saleem_reconstruction_attack | 6 | 0.7303 | 0.1952 | 0.0000 | 0.0000 | 0.0000 |
| 25 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | template_count_synthesis | 6 | 0.9515 | 0.0279 | 0.0000 | 0.0000 | 0.0000 |
| 25 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | independent_count | 6 | 0.2229 | 0.1209 | -0.0245 | 12.4007 | 0.9700 |
| 25 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | pairwise_refinement | 6 | 0.2229 | 0.1209 | -0.6752 | 12.4007 | 0.9700 |
| 25 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | saleem_reconstruction_attack | 6 | 0.2229 | 0.1209 | -0.5074 | 12.4007 | 0.9700 |
| 25 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | template_count_synthesis | 6 | 0.2229 | 0.1209 | -0.7285 | 12.4007 | 0.9700 |
| 25 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | independent_count | 6 | 0.2266 | 0.1253 | -0.0209 | 12.3087 | 0.9673 |
| 25 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | pairwise_refinement | 6 | 0.2284 | 0.1266 | -0.6698 | 12.3087 | 0.9673 |
| 25 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | saleem_reconstruction_attack | 6 | 0.2277 | 0.1251 | -0.5026 | 12.3087 | 0.9673 |
| 25 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | template_count_synthesis | 6 | 0.2275 | 0.1255 | -0.7239 | 12.3087 | 0.9673 |
| 25 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | independent_count | 6 | 0.2263 | 0.1208 | -0.0211 | 10.0103 | 0.9563 |
| 25 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | pairwise_refinement | 6 | 0.2516 | 0.1364 | -0.6465 | 10.0103 | 0.9563 |
| 25 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | saleem_reconstruction_attack | 6 | 0.2424 | 0.1205 | -0.4880 | 10.0103 | 0.9563 |
| 25 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | template_count_synthesis | 6 | 0.2463 | 0.1383 | -0.7051 | 10.0103 | 0.9563 |
| 50 | exact_count |  |  |  |  | independent_count | 18 | 0.2948 | 0.0825 |  | 0.0000 | 0.0000 |
| 50 | exact_count |  |  |  |  | pairwise_refinement | 6 | 0.8418 | 0.0285 |  | 0.0000 | 0.0000 |
| 50 | exact_count |  |  |  |  | saleem_reconstruction_attack | 6 | 0.6739 | 0.0664 |  | 0.0000 | 0.0000 |
| 50 | exact_count |  |  |  |  | template_count_synthesis | 18 | 0.8092 | 0.1078 |  | 0.0000 | 0.0000 |
| 50 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | independent_count | 9 | 0.2089 | 0.1158 | -0.0476 | 21.3309 | 0.9793 |
| 50 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | pairwise_refinement | 9 | 0.2280 | 0.1277 | -0.5989 | 21.3309 | 0.9793 |
| 50 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | saleem_reconstruction_attack | 9 | 0.2230 | 0.1198 | -0.4368 | 21.3309 | 0.9793 |
| 50 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | template_count_synthesis | 9 | 0.2295 | 0.1256 | -0.6343 | 21.3309 | 0.9793 |
| 50 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | independent_count | 18 | 0.2592 | 0.0851 | -0.0356 | 6.7342 | 0.9254 |
| 50 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | pairwise_refinement | 9 | 0.4297 | 0.1577 | -0.3908 | 6.8191 | 0.9438 |
| 50 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | saleem_reconstruction_attack | 9 | 0.3940 | 0.1218 | -0.2753 | 6.8191 | 0.9438 |
| 50 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | template_count_synthesis | 18 | 0.4495 | 0.1328 | -0.3597 | 6.7342 | 0.9254 |
| 50 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | independent_count | 18 | 0.2947 | 0.0828 | -0.0002 | 0.7758 | 0.5258 |
| 50 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | pairwise_refinement | 9 | 0.7203 | 0.0910 | -0.0636 | 0.7913 | 0.5344 |
| 50 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | saleem_reconstruction_attack | 9 | 0.6339 | 0.0665 | -0.0332 | 0.7913 | 0.5344 |
| 50 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | template_count_synthesis | 18 | 0.7579 | 0.0947 | -0.0513 | 0.7758 | 0.5258 |
| 50 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | independent_count | 18 | 0.2948 | 0.0825 | 0.0000 | 0.0000 | 0.0000 |
| 50 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | pairwise_refinement | 5 | 0.8454 | 0.0304 | 0.0000 | 0.0000 | 0.0000 |
| 50 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | saleem_reconstruction_attack | 5 | 0.6599 | 0.0634 | 0.0000 | 0.0000 | 0.0000 |
| 50 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | template_count_synthesis | 18 | 0.8092 | 0.1078 | 0.0000 | 0.0000 | 0.0000 |
| 50 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | independent_count | 18 | 0.1969 | 0.0857 | -0.0979 | 24.7519 | 0.9626 |
| 50 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | pairwise_refinement | 5 | 0.1865 | 0.1346 | -0.6588 | 24.7612 | 0.9916 |
| 50 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | saleem_reconstruction_attack | 5 | 0.1865 | 0.1346 | -0.4733 | 24.7612 | 0.9916 |
| 50 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | template_count_synthesis | 18 | 0.1969 | 0.0857 | -0.6123 | 24.7519 | 0.9626 |
| 50 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | independent_count | 18 | 0.2071 | 0.0864 | -0.0878 | 24.0226 | 0.9613 |
| 50 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | pairwise_refinement | 5 | 0.1920 | 0.1400 | -0.6534 | 24.0280 | 0.9936 |
| 50 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | saleem_reconstruction_attack | 5 | 0.1906 | 0.1368 | -0.4693 | 24.0280 | 0.9936 |
| 50 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | template_count_synthesis | 18 | 0.2084 | 0.0876 | -0.6008 | 24.0226 | 0.9613 |
| 50 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | independent_count | 18 | 0.2205 | 0.0872 | -0.0743 | 16.6959 | 0.9524 |
| 50 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | pairwise_refinement | 5 | 0.2347 | 0.1644 | -0.6106 | 16.5648 | 0.9836 |
| 50 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | saleem_reconstruction_attack | 5 | 0.2250 | 0.1498 | -0.4348 | 16.5648 | 0.9836 |
| 50 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | template_count_synthesis | 18 | 0.2512 | 0.1090 | -0.5581 | 16.6959 | 0.9524 |
| 100 | exact_count |  |  |  |  | independent_count | 3 | 0.1424 | 0.0113 |  | 0.0000 | 0.0000 |
| 100 | exact_count |  |  |  |  | pairwise_refinement | 3 | 0.7085 | 0.0183 |  | 0.0000 | 0.0000 |
| 100 | exact_count |  |  |  |  | saleem_reconstruction_attack | 3 | 0.7774 | 0.0314 |  | 0.0000 | 0.0000 |
| 100 | exact_count |  |  |  |  | template_count_synthesis | 3 | 0.9207 | 0.0131 |  | 0.0000 | 0.0000 |
| 100 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | independent_count | 9 | 0.2024 | 0.0974 | -0.0457 | 35.2071 | 0.9920 |
| 100 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | pairwise_refinement | 9 | 0.2468 | 0.1243 | -0.5909 | 35.2071 | 0.9920 |
| 100 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | saleem_reconstruction_attack | 9 | 0.2388 | 0.1063 | -0.6473 | 35.2071 | 0.9920 |
| 100 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | template_count_synthesis | 9 | 0.2523 | 0.1207 | -0.7885 | 35.2071 | 0.9920 |
| 100 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | independent_count | 9 | 0.2292 | 0.0840 | -0.0184 | 7.7767 | 0.9582 |
| 100 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | pairwise_refinement | 9 | 0.5169 | 0.1230 | -0.2996 | 7.7767 | 0.9582 |
| 100 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | saleem_reconstruction_attack | 9 | 0.5008 | 0.0973 | -0.3666 | 7.7767 | 0.9582 |
| 100 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | template_count_synthesis | 9 | 0.5718 | 0.1707 | -0.4962 | 7.7767 | 0.9582 |
| 100 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | independent_count | 9 | 0.2436 | 0.0802 | -0.0004 | 0.8240 | 0.5431 |
| 100 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | pairwise_refinement | 9 | 0.6670 | 0.0762 | -0.0216 | 0.8240 | 0.5431 |
| 100 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | saleem_reconstruction_attack | 9 | 0.6868 | 0.0765 | -0.0384 | 0.8240 | 0.5431 |
| 100 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | template_count_synthesis | 9 | 0.8331 | 0.1239 | -0.0725 | 0.8240 | 0.5431 |
| 100 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | independent_count | 3 | 0.1424 | 0.0113 | 0.0000 | 0.0000 | 0.0000 |
| 100 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | pairwise_refinement | 3 | 0.7085 | 0.0183 | 0.0000 | 0.0000 | 0.0000 |
| 100 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | saleem_reconstruction_attack | 3 | 0.7774 | 0.0314 | 0.0000 | 0.0000 | 0.0000 |
| 100 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | template_count_synthesis | 3 | 0.9207 | 0.0131 | 0.0000 | 0.0000 | 0.0000 |
| 100 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | independent_count | 3 | 0.0943 | 0.0222 | -0.0482 | 49.0667 | 0.9967 |
| 100 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | pairwise_refinement | 3 | 0.0940 | 0.0221 | -0.6145 | 49.0667 | 0.9967 |
| 100 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | saleem_reconstruction_attack | 3 | 0.0942 | 0.0223 | -0.6832 | 49.0667 | 0.9967 |
| 100 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | template_count_synthesis | 3 | 0.0943 | 0.0222 | -0.8264 | 49.0667 | 0.9967 |
| 100 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | independent_count | 3 | 0.0943 | 0.0252 | -0.0482 | 46.3533 | 0.9973 |
| 100 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | pairwise_refinement | 3 | 0.0974 | 0.0256 | -0.6111 | 46.3533 | 0.9973 |
| 100 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | saleem_reconstruction_attack | 3 | 0.0991 | 0.0273 | -0.6783 | 46.3533 | 0.9973 |
| 100 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | template_count_synthesis | 3 | 0.0977 | 0.0269 | -0.8229 | 46.3533 | 0.9973 |
| 100 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | independent_count | 3 | 0.1033 | 0.0196 | -0.0391 | 22.7367 | 0.9887 |
| 100 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | pairwise_refinement | 3 | 0.1543 | 0.0363 | -0.5542 | 22.7367 | 0.9887 |
| 100 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | saleem_reconstruction_attack | 3 | 0.1709 | 0.0449 | -0.6064 | 22.7367 | 0.9887 |
| 100 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | template_count_synthesis | 3 | 0.1363 | 0.0534 | -0.7844 | 22.7367 | 0.9887 |
| 250 | exact_count |  |  |  |  | independent_count | 3 | 0.1444 | 0.0042 |  | 0.0000 | 0.0000 |
| 250 | exact_count |  |  |  |  | pairwise_refinement | 3 | 0.5195 | 0.0049 |  | 0.0000 | 0.0000 |
| 250 | exact_count |  |  |  |  | saleem_reconstruction_attack | 3 | 0.6394 | 0.1040 |  | 0.0000 | 0.0000 |
| 250 | exact_count |  |  |  |  | template_count_synthesis | 3 | 0.9259 | 0.0054 |  | 0.0000 | 0.0000 |
| 250 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | independent_count | 9 | 0.2191 | 0.0850 | -0.0263 | 56.3758 | 0.9953 |
| 250 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | pairwise_refinement | 9 | 0.3035 | 0.1189 | -0.3420 | 56.3758 | 0.9953 |
| 250 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | saleem_reconstruction_attack | 9 | 0.3139 | 0.1025 | -0.4364 | 56.3758 | 0.9953 |
| 250 | dp_per_query_eps0.01 | per_query | 0.01 | 0.01 | 5 | template_count_synthesis | 9 | 0.3543 | 0.1313 | -0.6978 | 56.3758 | 0.9953 |
| 250 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | independent_count | 9 | 0.2441 | 0.0780 | -0.0030 | 9.1058 | 0.9613 |
| 250 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | pairwise_refinement | 9 | 0.5111 | 0.0615 | -0.0627 | 9.1058 | 0.9613 |
| 250 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | saleem_reconstruction_attack | 9 | 0.5995 | 0.0928 | -0.1228 | 9.1058 | 0.9613 |
| 250 | dp_per_query_eps0.1 | per_query | 0.1 | 0.1 | 50 | template_count_synthesis | 9 | 0.7327 | 0.1259 | -0.2671 | 9.1058 | 0.9613 |
| 250 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | independent_count | 9 | 0.2482 | 0.0779 | 0.0008 | 0.8324 | 0.5447 |
| 250 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | pairwise_refinement | 9 | 0.5378 | 0.0469 | 0.0022 | 0.8324 | 0.5447 |
| 250 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | saleem_reconstruction_attack | 9 | 0.6738 | 0.0678 | 0.0005 | 0.8324 | 0.5447 |
| 250 | dp_per_query_eps1 | per_query | 1 | 1 | 500 | template_count_synthesis | 9 | 0.8681 | 0.1107 | -0.0203 | 0.8324 | 0.5447 |
| 250 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | independent_count | 3 | 0.1444 | 0.0042 | 0.0000 | 0.0000 | 0.0000 |
| 250 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | pairwise_refinement | 3 | 0.5195 | 0.0049 | 0.0000 | 0.0000 | 0.0000 |
| 250 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | saleem_reconstruction_attack | 3 | 0.6394 | 0.1040 | 0.0000 | 0.0000 | 0.0000 |
| 250 | dp_per_query_eps10 | per_query | 10 | 10 | 5000 | template_count_synthesis | 3 | 0.9259 | 0.0054 | 0.0000 | 0.0000 | 0.0000 |
| 250 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | independent_count | 3 | 0.1063 | 0.0179 | -0.0382 | 121.1907 | 1.0000 |
| 250 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | pairwise_refinement | 3 | 0.1059 | 0.0178 | -0.4135 | 121.1907 | 1.0000 |
| 250 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | saleem_reconstruction_attack | 3 | 0.1065 | 0.0182 | -0.5329 | 121.1907 | 1.0000 |
| 250 | dp_composed_eps0.1 | composed | 0.1 | 0.0002 | 0.1 | template_count_synthesis | 3 | 0.1067 | 0.0180 | -0.8192 | 121.1907 | 1.0000 |
| 250 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | independent_count | 3 | 0.1083 | 0.0210 | -0.0361 | 102.0547 | 0.9993 |
| 250 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | pairwise_refinement | 3 | 0.1147 | 0.0199 | -0.4048 | 102.0547 | 0.9993 |
| 250 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | saleem_reconstruction_attack | 3 | 0.1202 | 0.0282 | -0.5192 | 102.0547 | 0.9993 |
| 250 | dp_composed_eps1 | composed | 1 | 0.002 | 1 | template_count_synthesis | 3 | 0.1237 | 0.0279 | -0.8022 | 102.0547 | 0.9993 |
| 250 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | independent_count | 3 | 0.1246 | 0.0097 | -0.0199 | 29.6387 | 0.9907 |
| 250 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | pairwise_refinement | 3 | 0.2575 | 0.0150 | -0.2620 | 29.6387 | 0.9907 |
| 250 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | saleem_reconstruction_attack | 3 | 0.2847 | 0.0517 | -0.3547 | 29.6387 | 0.9907 |
| 250 | dp_composed_eps10 | composed | 10 | 0.02 | 10 | template_count_synthesis | 3 | 0.2934 | 0.0359 | -0.6324 | 29.6387 | 0.9907 |

## Claim Boundary

Per-query epsilon is the privacy loss of one count. It is not the privacy loss of the 500-SNP release.
Composed epsilon is the basic-composition budget of the whole count vector.
Rows that match the noisy counts can still differ from the true genotypes; F1 measures that difference.
