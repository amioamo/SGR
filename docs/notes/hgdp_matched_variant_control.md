# Matched-Variant Within-Cohort Control

The variant set in each window is the set retained for the HGDP European target.
Both arms use the same 1000 Genomes EUR public reference. The within-cohort arm
replaces the HGDP beacon with the prototype 1000 Genomes beacon.

## Configuration

- beacon sizes: 25, 50, 100
- reference panel size: 100
- replicates: 3
- minimum kept variants: 100

## Carrier F1 on the Shared Sites

| N | method | target | n | mean F1 | median F1 | mean Hamming | Sim ≥ 0.95 | invalid |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 25 | reference_only | 1kg_target | 15 | 0.782939 | 0.769481 | 489.60 | 0.626667 | 15 |
| 25 | reference_only | hgdp_target | 15 | 0.727207 | 0.710630 | 601.53 | 0.640000 | 15 |
| 25 | template_count_synthesis | 1kg_target | 15 | 0.863779 | 0.823697 | 252.40 | 0.728000 | 0 |
| 25 | template_count_synthesis | hgdp_target | 15 | 0.850069 | 0.796610 | 266.53 | 0.722667 | 0 |
| 50 | reference_only | 1kg_target | 15 | 0.769634 | 0.777362 | 885.47 | 0.710667 | 15 |
| 50 | reference_only | hgdp_target | 15 | 0.762146 | 0.789969 | 863.47 | 0.725333 | 15 |
| 50 | template_count_synthesis | 1kg_target | 15 | 0.877274 | 0.858079 | 451.07 | 0.760000 | 0 |
| 50 | template_count_synthesis | hgdp_target | 15 | 0.877857 | 0.864635 | 428.00 | 0.790667 | 0 |
| 100 | reference_only | 1kg_target | 15 | 0.823731 | 0.847702 | 1358.93 | 0.790667 | 15 |
| 100 | reference_only | hgdp_target | 15 | 0.834883 | 0.837637 | 1316.67 | 0.776667 | 15 |
| 100 | template_count_synthesis | 1kg_target | 15 | 0.892628 | 0.864855 | 805.60 | 0.812667 | 0 |
| 100 | template_count_synthesis | hgdp_target | 15 | 0.888566 | 0.863506 | 840.40 | 0.794667 | 0 |

## Paired Difference, HGDP Target minus 1KG Target

Positive means the HGDP beacon is reconstructed more accurately than the 1KG beacon
on the same sites and the same reference panel.

| N | method | metric | n | mean within | mean cross | mean delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |
|---:|---|---|---:|---:|---:|---:|---|---:|---:|---:|
| 25 | reference_only | f1 | 15 | 0.782939 | 0.727207 | -0.055732 | [-0.112388, 0.003211] | 3/0/12 | 0.0352 | 0.0637 |
| 25 | reference_only | proportion_similarity_ge_095 | 15 | 0.626667 | 0.640000 | 0.013333 | [-0.042667, 0.064000] | 8/1/6 | 0.7905 | 0.4126 |
| 25 | template_count_synthesis | f1 | 15 | 0.863779 | 0.850069 | -0.013710 | [-0.040891, 0.010585] | 6/0/9 | 0.6072 | 0.3894 |
| 25 | template_count_synthesis | proportion_similarity_ge_095 | 15 | 0.728000 | 0.722667 | -0.005333 | [-0.064000, 0.050667] | 5/6/4 | 1.0000 | 0.9336 |
| 50 | reference_only | f1 | 15 | 0.769634 | 0.762146 | -0.007489 | [-0.071095, 0.051008] | 8/0/7 | 1.0000 | 0.9341 |
| 50 | reference_only | proportion_similarity_ge_095 | 15 | 0.710667 | 0.725333 | 0.014667 | [-0.022667, 0.054667] | 9/2/4 | 0.2668 | 0.6221 |
| 50 | template_count_synthesis | f1 | 15 | 0.877274 | 0.877857 | 0.000583 | [-0.016065, 0.020707] | 7/0/8 | 1.0000 | 0.9341 |
| 50 | template_count_synthesis | proportion_similarity_ge_095 | 15 | 0.760000 | 0.790667 | 0.030667 | [-0.033333, 0.102667] | 6/6/3 | 0.5078 | 0.5859 |
| 100 | reference_only | f1 | 15 | 0.823731 | 0.834883 | 0.011152 | [-0.033476, 0.078180] | 7/0/8 | 1.0000 | 0.5614 |
| 100 | reference_only | proportion_similarity_ge_095 | 15 | 0.790667 | 0.776667 | -0.014000 | [-0.032667, 0.004667] | 4/0/11 | 0.1185 | 0.1108 |
| 100 | template_count_synthesis | f1 | 15 | 0.892628 | 0.888566 | -0.004062 | [-0.012551, 0.004820] | 6/0/9 | 0.6072 | 0.2769 |
| 100 | template_count_synthesis | proportion_similarity_ge_095 | 15 | 0.812667 | 0.794667 | -0.018000 | [-0.047333, 0.010667] | 3/4/8 | 0.2266 | 0.2764 |

## Full Method by Window

| window | sites | N | within F1 | cross F1 | cross − within |
|---|---:|---:|---:|---:|---:|
| chr1_12795449_12810845_dense500_w0 | 298 | 25 | 0.987876 | 0.989597 | 0.001721 |
| chr1_12795449_12810845_dense500_w0 | 298 | 50 | 0.982016 | 0.990036 | 0.008020 |
| chr1_12795449_12810845_dense500_w0 | 298 | 100 | 0.984778 | 0.988335 | 0.003557 |
| chr1_246210538_246267706_dense500_w1 | 408 | 25 | 0.979868 | 0.979902 | 0.000034 |
| chr1_246210538_246267706_dense500_w1 | 408 | 50 | 0.986154 | 0.981424 | -0.004730 |
| chr1_246210538_246267706_dense500_w1 | 408 | 100 | 0.985112 | 0.975315 | -0.009798 |
| chr20_30915197_30991425_dense500_w1 | 211 | 25 | 0.830833 | 0.784496 | -0.046337 |
| chr20_30915197_30991425_dense500_w1 | 211 | 50 | 0.831550 | 0.851052 | 0.019502 |
| chr20_30915197_30991425_dense500_w1 | 211 | 100 | 0.855111 | 0.865755 | 0.010644 |
| chr22_11560143_11630262_dense500_w0 | 405 | 25 | 0.707028 | 0.739948 | 0.032921 |
| chr22_11560143_11630262_dense500_w0 | 405 | 50 | 0.761673 | 0.755455 | -0.006217 |
| chr22_11560143_11630262_dense500_w0 | 405 | 100 | 0.789239 | 0.778533 | -0.010706 |
| chr22_22346872_22422653_dense500_w1 | 204 | 25 | 0.813289 | 0.756401 | -0.056887 |
| chr22_22346872_22422653_dense500_w1 | 204 | 50 | 0.824977 | 0.811318 | -0.013660 |
| chr22_22346872_22422653_dense500_w1 | 204 | 100 | 0.848901 | 0.834894 | -0.014008 |

## Claim Boundary

This comparison fixes the harmonized site set. It does not restore the original
500-SNP windows. The two beacons are different people, so the paired difference
mixes cohort shift with the particular beacon draw.

