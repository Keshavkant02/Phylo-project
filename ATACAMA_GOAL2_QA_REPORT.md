# Atacama Goal 2 QA Report

## What changed

- Rebuilt `soil_microbiome_16s_class_safe_colab.ipynb` as an Atacama-only 16S ASV notebook.
- Removed the synthetic `Soil_ASV_A` / `Soil_ASV_B` teaching path from the student-facing notebook.
- Removed visible form controls, live BLAST, neighbor joining, IQ-TREE, and bootstrap content from the student-facing notebook.
- Added `build_atacama_soil_asv_phylogeny_colab.py` and `verify_atacama_soil_asv_phylogeny_colab.py`.

## Data source

The generated cache uses the local QIIME 2 Atacama artifacts in `tmp/atacama_qiime2_source`:

- `atacama-table.qza`
- `atacama-rep-seqs.qza`
- `sample_metadata.tsv`

The generated table has 61 samples and 401 ASVs.

## Scientific honesty notes

- No valid SILVA taxonomy artifact was present in the local source directory. The notebook therefore displays `Unassigned at genus level` instead of inventing taxonomy.
- The 61-sample source table has only nine ASVs present in at least 10% of samples. The q-value section still tests the requested top 50 ASVs by prevalence, but the lollipop plot labels only significant ASVs that meet the 10% prevalence threshold.
- BH-adjusted q-values are computed live from the loaded abundance table and metadata using CLR-transformed abundance, a 0.5 pseudo-count, and `statsmodels.stats.multitest.multipletests(method="fdr_bh")`.

## Verification

Ran:

```text
python build_atacama_soil_asv_phylogeny_colab.py
python verify_atacama_soil_asv_phylogeny_colab.py
```

Verifier results:

- 35 total notebook cells.
- 14 code cells and 21 markdown cells.
- 0 code-cell execution errors.
- 0 visible Colab form controls.
- 0 synthetic teaching sections in the notebook source.
- 61 samples loaded.
- 50 ASVs tested for q-values.
- 20 ASVs used for relative abundance plots.
- 12 ASVs used for the UPGMA tree.
- 8 ASVs used for the alignment heatmap.
- 9 figure outputs saved under `soil_16s_class_cache/goal2_figure_checks`.
- 4 styled table outputs generated.

Visual spot checks passed for:

- Alignment heatmap conserved-column desaturation.
- UPGMA tree label readability.
- Distance matrix palette and labels.
- Relative abundance plot layout.
- Alpha diversity scatter plots.
- Raw p-value versus q-value plot.
- Lollipop plot sort/readability after prevalence filtering.
