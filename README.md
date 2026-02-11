# Single-Cell Reanalysis

## Smoke test (end-to-end)
This repo includes a minimal, end-to-end smoke test to confirm the pipeline runs start -> finish.

### 1) Create environment (Windows / conda)
    conda create -n scr_smoke -c conda-forge python=3.11 scanpy anndata python-igraph leidenalg pyyaml matplotlib pandas numpy -y
    conda activate scr_smoke

### 2) Generate toy data
    python src/pipeline/make_toy_data.py

### 3) Run pipeline
    python src/pipeline/run.py --config configs/smoke.yaml

### Outputs
- outputs/processed_smoke.h5ad
- reports/smoke/qc_summary.csv
- reports/smoke/qc_violin.png
- reports/smoke/qc_counts_vs_mito.png
- reports/smoke/umap_leiden.png

## Full analysis run (discovery + validation figures)
From the repo root, with the analysis environment activated:

    conda activate scr_smoke
    python scripts/run_full.py --config configs/full.yaml

This regenerates:
- results/tables/state_summary.tsv
- results/tables/paired_subject_deltas.tsv
- results/figures/fig1*.png
- results/figures/figV*.png (if validation tables exist)
