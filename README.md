# SGR: Structure-Guided Reconstruction Attacks on Genomic Data Sharing Beacons

Reference implementation of a structure-guided attack (SGR) that reconstructs anonymous participant-level carrier matrices from aggregate Beacon carrier-count responses using public reference haplotypes.

The attack entry point is `reference_haplotype_reconstruction` in `src/beacon_reconstruction/reconstruction.py`. It greedily selects reference-derived multi-locus templates while tracking the released carrier-frequency profile, permits template reuse, and then projects the selected template onto the released count constraints exactly.

Associated manuscript: *Structure-Guided Reconstruction Attacks on Genomic Data Sharing Beacons*.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Python 3.11+ recommended. PyTorch is required only for the LD/AF optimization baseline (Saleem et al.).

## Data

Preprocessed 1000 Genomes panels and manifests live under `data/`. See [data/README.md](data/README.md) for file layout, sample partitions, and instructions to rebuild panels from IGSR downloads.

## Experiments

Run commands from the repository root. Outputs are CSV tables and a markdown note under `results/`. See [docs/experiments.md](docs/experiments.md).

```bash
python3 scripts/run_cohort_scaling.py --help
python3 scripts/run_release_policies.py --help
python3 scripts/run_participant_privacy.py --help
python3 scripts/run_runtime_benchmark.py --help
python3 scripts/run_reference_matching.py --help
python3 scripts/run_multiancestry_reference_matching.py --help
python3 scripts/run_template_stage_ablation.py --help
python3 scripts/run_greedy_template_optimality.py --help
python3 scripts/run_nested_reference_panel_size.py --help
python3 scripts/run_ld_span_query_evaluation.py --help
python3 scripts/run_dp_count_release.py --help
python3 scripts/run_partial_genotype_linkage.py --help
python3 scripts/run_hgdp_eur_cross_dataset.py --help
python3 scripts/run_hgdp_matched_variant_control.py --help
python3 scripts/run_hgdp_eur_paper_methods.py --help
```

## Core library

Baselines and evaluation utilities are in the same package:

- `src/beacon_reconstruction/baselines.py` — count-only and local LD refinement
- `src/beacon_reconstruction/ld_af_baseline.py` — LD/AF optimization baseline
- `src/beacon_reconstruction/evaluation.py` — row alignment and privacy metrics

## Citation

If you use this code, please cite the associated article and archive DOI (to be added upon Zenodo deposit).

## License

MIT License. See [LICENSE](LICENSE).
