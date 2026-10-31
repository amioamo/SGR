# Experiments

Run these from the repository root after installing `requirements.txt`. Each script writes CSV tables and a markdown note under `results/`. Completed-run notes are in [notes/](notes/).

SGR is `reference_haplotype_reconstruction`: pick reference rows by a proportional greedy rule that tracks the released carrier counts, allowing row reuse, then project that template onto the counts exactly.

Shared defaults for the 1000 Genomes runs are the six 500-SNP windows, three replicates, and the shared reference panel of 100 individuals disjoint from the beacon.

## Stage ablation

`scripts/run_template_stage_ablation.py`

Compares the full attack with three removals of its stages: reference rows sampled without using the counts (`reference_only`), those rows after exact projection (`random_template_projection`), and the greedy template without projection (`greedy_template_without_projection`).

```bash
python3 scripts/run_template_stage_ablation.py \
  --beacon-sizes 25,50,100,250 \
  --max-windows 6 \
  --max-replicates 3 \
  --shared-reference-panel
```

Default output: `results/block_inheritance_v1/template_stage_ablation_shared_panel/`.

On the completed run, exact projection is what makes the reconstruction count-feasible. Random templates stay well below the greedy template after projection. Skipping projection leaves a small carrier-count error and is marked invalid for that reason. Note: [notes/template_stage_ablation.md](notes/template_stage_ablation.md).

## Greedy template versus the exact L1 template

`scripts/run_greedy_template_optimality.py`

For small reference size `R` and beacon size `N`, enumerates every multiset of `N` reference rows and compares the greedy Stage-1 template with the template that minimizes `||y - c(T)||_1`.

```bash
python3 scripts/run_greedy_template_optimality.py
```

Default settings are `R:N` in `5:5,6:4,8:5,8:8,10:6,10:8`. Default output: `results/block_inheritance_v1/greedy_template_optimality/`.

The completed run finds a small gap. At `R=5, N=5` the mean absolute L1 gap is 1.34 and the optimal-hit rate is 0.90. At `R=10, N=8` (24,310 templates) the mean absolute gap is 2.42, the mean relative gap is 0.037, and the hit rate is 0.62. Note: [notes/greedy_template_optimality.md](notes/greedy_template_optimality.md).

## Nested reference-panel size

`scripts/run_nested_reference_panel_size.py`

Holds the beacon fixed and grows the EUR reference as a nested prefix, so each larger panel contains the smaller ones. Panel sizes are 10, 25, 50, 100, 200, and 400.

```bash
python3 scripts/run_nested_reference_panel_size.py
```

Default output: `results/block_inheritance_v1/reference_panel_size_nested_eur/`.

At beacon size 50, mean F1 rises from 0.705 at `R=10` to 0.815 at `R=100` and 0.821 at `R=400`. Direct sample overlap is 0. Note: [notes/nested_reference_panel_size.md](notes/nested_reference_panel_size.md).

## LD, density, and genomic span

`scripts/run_ld_span_query_evaluation.py`

Uses the 36 prepared panels in `data/ld_span_queries/`. Every panel has 200 SNPs. Regimes are dense, LD-pruned at r^2 0.80 / 0.50 / 0.10, sparse (50 kb gap), and long-range (about 10 Mb). These regions sit outside the six main windows.

```bash
python3 scripts/run_ld_span_query_evaluation.py \
  --manifest data/ld_span_queries/ld_span_query_manifest.csv \
  --beacon-sizes 25,50 \
  --max-replicates 3 \
  --no-run-saleem
```

Add the LD/AF baseline with `--saleem-only --merge-metrics <metrics.csv>` on a second pass, or use `scripts/run_ld_span_query_experiment.sh` when `plink2` and the upstream candidate pfiles are available. Rebuilding panels: `python3 scripts/prepare_ld_span_query_panels.py --plink2 /path/to/plink2`.

Default output: `results/block_inheritance_v1/ld_span_queries_shared_panel/`.

At `N=50`, matched-template mean F1 is 0.828 on dense panels, 0.571 / 0.515 / 0.437 at the three r^2 caps, and about 0.42 on sparse and long-range panels. The template advantage over local LD and over LD/AF optimization is concentrated in the short dense panels. Note: [notes/ld_span_queries.md](notes/ld_span_queries.md).

## Discrete Laplace count release

`scripts/run_dp_count_release.py`

Adds independent discrete Laplace noise of sensitivity 1 to each carrier count and clips the result to `[0, N]`. Two budgets are reported. Per-query accounting spends epsilon on every SNP, so `m` SNPs compose to `m * epsilon`. Composed accounting splits epsilon across the SNPs, so the joint release has budget epsilon.

```bash
python3 scripts/run_dp_count_release.py \
  --beacon-sizes 25,50,100,250 \
  --epsilons 0.01,0.1,1,10
```

Default output: `results/block_inheritance_v1/dp_count_release_shared_panel/`.

Exact counts preserve SGR. A joint budget of epsilon 0.1 or 1, split across 500 SNPs, drops SGR mean F1 at `N=50` from 0.809 to about 0.20. Spending epsilon 1 on every SNP leaves more of the signal (mean F1 0.758 at `N=50`). Note: [notes/dp_count_release.md](notes/dp_count_release.md).

## Partial-genotype linkage

`scripts/run_partial_genotype_linkage.py`

After row alignment, treats a subset of the true carrier calls as known and asks how often the matching reconstructed row is the nearest, or among the five nearest, reconstructed rows.

```bash
python3 scripts/run_partial_genotype_linkage.py \
  --beacon-sizes 25,50,100,250 \
  --known-snp-counts 10,25,50,100
```

Default output: `results/block_inheritance_v1/partial_genotype_linkage_shared_panel/`.

At `N=25` and 10 known SNPs, the matched-template top-1 linkage rate is 0.167 and top-5 is 0.518, against 0.044 and 0.244 for independent counts. Both rise as more SNPs are known. Note: [notes/partial_genotype_linkage.md](notes/partial_genotype_linkage.md).

## HGDP European target

Three scripts use `data/external/hgdp_1kg/`. The beacon is the HGDP European subset of the gnomAD HGDP+1KG callset, harmonized to the same target allele as the local 1000 Genomes window. The reference is 1000 Genomes EUR and does not overlap the beacon. Window `chr20:29757801-29829649` has no overlapping HGDP variants and is skipped, leaving five windows.

`scripts/prepare_hgdp_eur_target.py` rebuilds those slices from the public gnomAD VCF. It needs `bcftools`.

Cross-dataset reconstruction:

```bash
python3 scripts/run_hgdp_eur_cross_dataset.py
```

Within-cohort control on the same harmonized sites, so the target is 1000 Genomes rather than HGDP:

```bash
python3 scripts/run_hgdp_matched_variant_control.py
```

Count, LD/AF optimization, and SGR on the HGDP beacon:

```bash
python3 scripts/run_hgdp_eur_paper_methods.py
```

Outputs:

- `results/block_inheritance_v1/hgdp_eur_cross_dataset/`
- `results/block_inheritance_v1/hgdp_matched_variant_control/`
- `results/block_inheritance_v1/hgdp_eur_paper_methods/`

