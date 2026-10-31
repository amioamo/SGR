# HGDP European Target, 1000 Genomes EUR Reference

The beacon is a sample of HGDP Europeans from the gnomAD v3.1.2 HGDP+1KG
callset. The reference is the matched 1000 Genomes EUR public panel used in
the main experiment. Sites are restricted to variants present in both, with the
1000 Genomes target allele, and called in every retained HGDP sample.
Samples below a 95% call rate on the matched sites are removed before that filter.

## Configuration

- HGDP European samples: 146
- beacon sizes: 25, 50, 100
- reference panel size: 100
- windows requested: 6
- windows analyzed: 5
- replicates: 3

## Variant Overlap

| window | 1KG variants | kept | exact | flipped | missing site | dropped missing GT |
|---|---:|---:|---:|---:|---:|---:|
| chr1_12795449_12810845_dense500_w0 | 500 | 298 | 298 | 0 | 2 | 200 |
| chr1_246210538_246267706_dense500_w1 | 500 | 408 | 408 | 0 | 2 | 90 |
| chr20_29757801_29829649_dense500_w0 | 500 | 0 | 0 | 0 | 500 | 0 |
| chr20_30915197_30991425_dense500_w1 | 500 | 211 | 211 | 0 | 4 | 285 |
| chr22_11560143_11630262_dense500_w0 | 500 | 405 | 405 | 0 | 5 | 90 |
| chr22_22346872_22422653_dense500_w1 | 500 | 204 | 204 | 0 | 5 | 291 |

## Carrier F1

| N | method | n | mean F1 | median F1 | mean Hamming | mean Sim | Sim ≥ 0.95 | invalid |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 25 | reference_only | 15 | 0.727207 | 0.710630 | 601.53 | 0.918945 | 0.640000 | 15 |
| 25 | template_count_synthesis | 15 | 0.850069 | 0.796610 | 266.53 | 0.959519 | 0.722667 | 0 |
| 50 | reference_only | 15 | 0.762146 | 0.789969 | 863.47 | 0.940940 | 0.725333 | 15 |
| 50 | template_count_synthesis | 15 | 0.877857 | 0.864635 | 428.00 | 0.968120 | 0.790667 | 0 |
| 100 | reference_only | 15 | 0.834883 | 0.837637 | 1316.67 | 0.955026 | 0.776667 | 15 |
| 100 | template_count_synthesis | 15 | 0.888566 | 0.863506 | 840.40 | 0.969359 | 0.794667 | 0 |

## Paired Difference, Full Method minus Reference-only

| N | metric | n | mean delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |
|---:|---|---:|---:|---|---:|---:|---:|
| 25 | f1 | 15 | 0.122862 | [0.073851, 0.180321] | 15/0/0 | 6.10e-05 | 6.10e-05 |
| 25 | proportion_similarity_ge_095 | 15 | 0.082667 | [0.034667, 0.133333] | 10/3/2 | 0.0386 | 0.0049 |
| 50 | f1 | 15 | 0.115711 | [0.053699, 0.205498] | 13/0/2 | 0.0074 | 0.0006 |
| 50 | proportion_similarity_ge_095 | 15 | 0.065333 | [0.013333, 0.116000] | 12/0/3 | 0.0352 | 0.0170 |
| 100 | f1 | 15 | 0.053684 | [0.030932, 0.078066] | 14/0/1 | 0.0010 | 0.0006 |
| 100 | proportion_similarity_ge_095 | 15 | 0.018000 | [-0.010667, 0.046667] | 10/1/4 | 0.1796 | 0.1867 |

## Claim Boundary

HGDP and 1000 Genomes are different sample collections. N=250 is not run because
the retained HGDP European subset has 146 individuals. Windows with fewer than
100 harmonized variants are omitted; chr20:29757801-29829649 has no variants in
this gnomAD VCF. Reference-only draws rows from the 1000 Genomes panel and does
not use the beacon counts.

