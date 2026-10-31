# HGDP European Target: Count, LD/AF, and SGR

The beacon is the retained HGDP European cohort. The reference is the 1000 Genomes
EUR public panel. Sites are the harmonized variants used in the cross-dataset run.
Count-based assignment uses only the released counts. LD/AF optimization is the
Saleem attack with its published initialization. SGR is template-count synthesis.

## Configuration

- HGDP samples: 146
- beacon sizes: 25, 50, 100
- windows analyzed: 5
- replicates: 3

## Carrier F1 and High-Similarity Coverage

| N | method | n | mean F1 | SD | Sim ≥ 0.95 | mean Hamming | invalid |
|---:|---|---:|---:|---:|---:|---:|---:|
| 25 | Count-based | 15 | 0.318606 | 0.097199 | 0.002667 | 1551.07 | 0 |
| 25 | LD/AF optimization | 15 | 0.760516 | 0.106002 | 0.477333 | 555.47 | 0 |
| 25 | SGR | 15 | 0.850069 | 0.116826 | 0.722667 | 266.53 | 0 |
| 50 | Count-based | 15 | 0.309687 | 0.098943 | 0.004000 | 2922.67 | 0 |
| 50 | LD/AF optimization | 15 | 0.751398 | 0.118422 | 0.486667 | 1180.93 | 0 |
| 50 | SGR | 15 | 0.877857 | 0.098296 | 0.790667 | 428.00 | 0 |
| 100 | Count-based | 15 | 0.308284 | 0.098823 | 0.004667 | 6122.80 | 0 |
| 100 | LD/AF optimization | 15 | 0.738551 | 0.098622 | 0.486000 | 2544.00 | 0 |
| 100 | SGR | 15 | 0.888566 | 0.084988 | 0.794667 | 840.40 | 0 |

## Paired Difference, SGR minus Baseline

| N | baseline | metric | n | mean SGR | mean baseline | mean delta | 95% CI | wins/ties/losses | sign p | wilcoxon p |
|---:|---|---|---:|---:|---:|---:|---|---:|---:|---:|
| 25 | Count-based | f1 | 15 | 0.850069 | 0.318606 | 0.531463 | [0.442895, 0.629035] | 15/0/0 | 6.10e-05 | 6.10e-05 |
| 25 | Count-based | proportion_similarity_ge_095 | 15 | 0.722667 | 0.002667 | 0.720000 | [0.597333, 0.840000] | 15/0/0 | 6.10e-05 | 0.0006 |
| 25 | LD/AF optimization | f1 | 15 | 0.850069 | 0.760516 | 0.089553 | [0.054253, 0.129369] | 15/0/0 | 6.10e-05 | 6.10e-05 |
| 25 | LD/AF optimization | proportion_similarity_ge_095 | 15 | 0.722667 | 0.477333 | 0.245333 | [0.157333, 0.341400] | 14/1/0 | 0.0001 | 0.0010 |
| 50 | Count-based | f1 | 15 | 0.877857 | 0.309687 | 0.568170 | [0.481154, 0.658147] | 15/0/0 | 6.10e-05 | 6.10e-05 |
| 50 | Count-based | proportion_similarity_ge_095 | 15 | 0.790667 | 0.004000 | 0.786667 | [0.670667, 0.886667] | 15/0/0 | 6.10e-05 | 0.0006 |
| 50 | LD/AF optimization | f1 | 15 | 0.877857 | 0.751398 | 0.126459 | [0.074924, 0.184096] | 14/0/1 | 0.0010 | 0.0001 |
| 50 | LD/AF optimization | proportion_similarity_ge_095 | 15 | 0.790667 | 0.486667 | 0.304000 | [0.169333, 0.450667] | 13/1/1 | 0.0018 | 0.0012 |
| 100 | Count-based | f1 | 15 | 0.888566 | 0.308284 | 0.580282 | [0.504123, 0.662850] | 15/0/0 | 6.10e-05 | 6.10e-05 |
| 100 | Count-based | proportion_similarity_ge_095 | 15 | 0.794667 | 0.004667 | 0.790000 | [0.690000, 0.882000] | 15/0/0 | 6.10e-05 | 0.0007 |
| 100 | LD/AF optimization | f1 | 15 | 0.888566 | 0.738551 | 0.150015 | [0.095677, 0.211638] | 15/0/0 | 6.10e-05 | 6.10e-05 |
| 100 | LD/AF optimization | proportion_similarity_ge_095 | 15 | 0.794667 | 0.486000 | 0.308667 | [0.194000, 0.435333] | 14/1/0 | 0.0001 | 0.0001 |

## Claim Boundary

HGDP and 1000 Genomes participants do not overlap. Both cohorts are harmonized
through the gnomAD v3.1.2 HGDP+1KG callset at the retained sites, so this is a
cross-cohort comparison rather than an independent variant-calling comparison.
Windows with fewer than 100 harmonized variants are omitted.

