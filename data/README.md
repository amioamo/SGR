# Data for reference haplotype reconstruction experiments

This directory contains manifests and preprocessed genotype panels used by the manuscript experiments.

## Source data

Genotypes are derived from the [1000 Genomes Project](https://www.internationalgenome.org/) (Phase 3, GRCh38). Download VCF/PLINK files from the International Genome Sample Resource and place processed carrier matrices in the paths listed below.

## Directory layout

```
data/
├── manifests/
│   ├── 1000g_eur_chr22_prototype_samples.csv   # sample IDs and replicate roles
│   ├── 1000g_eur_chr22_prototype_seeds.csv     # three fixed partition replicates
│   └── 1000g_eur_chr22_prototype_variants.csv  # prototype variant metadata
└── multiwindow/
    ├── multiwindow_manifest.csv                # six 500-SNP query windows
    ├── plink/                                  # EUR carrier matrices (.raw)
    ├── wrong_population/                       # AFR control matrices for mechanism tests
    ├── ancestry_panels/
    │   ├── eur/                                # matched EUR reference rows
    │   ├── afr/
    │   ├── eas/
    │   ├── sas/
    │   └── amr/                                # mismatched-ancestry reference rows
    └── snplists/                               # SNP lists per window (optional)
```

## Query windows

Six consecutive 500-SNP windows on chromosomes 1, 20, and 22 (two windows per chromosome, ≥1 Mb apart). See `multiwindow/multiwindow_manifest.csv` for coordinates and file paths.

Each `.raw` file is a PLINK `.raw` matrix converted to binary carrier states: for each variant, the less frequent allele in EUR filtering genotypes is the target allele; an individual is coded as carrier if at least one copy is present.

## Sample partitions

For each of three fixed replicates (`seeds.csv`):

- **Beacon cohort**: nested sizes N ∈ {25, 50, 100, 250}
- **Shared reference panel**: R = 100 disjoint individuals (`--shared-reference-panel`)
- Beacon and reference individuals do not overlap within a replicate

Population-matching experiments mix EUR and one mismatched super-population (AFR, EAS, SAS, or AMR) in the auxiliary reference panel at fractions 0%, 25%, 50%, 75%, and 100%.

## Building panels from 1000 Genomes

If `.raw` files are not bundled with the release, rebuild them with:

```bash
python3 scripts/prep/prepare_multiwindow_panels.py --help
```

Then copy or symlink the generated `.raw` files into `data/multiwindow/plink/` (EUR beacon windows), `data/multiwindow/wrong_population/` (AFR control windows), and `data/multiwindow/ancestry_panels/{eur,afr,eas,sas,amr}/` (reference panels), matching the paths in `multiwindow_manifest.csv`.

This release bundle includes the preprocessed `.raw` panels (~22 MB total): 6 EUR beacon windows, 6 AFR wrong-population controls, and 6 windows × 5 ancestry reference panels.

## Additional panels

`ld_span_queries/` holds 36 prepared 200-SNP panels on chromosomes 1, 20, and 22. Each region has a dense panel, three LD-pruned panels (`r^2` caps 0.80, 0.50, and 0.10), a sparse panel (50 kb minimum gap), and a long-range panel (about 10 Mb). See `ld_span_queries/ld_span_query_manifest.csv`. Rebuilding these panels from candidate pfiles requires `plink2` (`scripts/prepare_ld_span_query_panels.py --plink2`).

`external/hgdp_1kg/` holds HGDP European genotypes harmonized onto the six main windows, plus sample lists. One window, `chr20:29757801-29829649`, has no overlapping HGDP variants. Rebuild instructions are in `scripts/prepare_hgdp_eur_target.py`.

## Outputs

Experiment scripts write CSV metrics and a markdown note under `results/`. See [docs/experiments.md](../docs/experiments.md) for each script and where it writes.
