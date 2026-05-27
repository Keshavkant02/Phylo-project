# Atacama Goal 2 QA Report

## What changed

- Rebuilt `soil_microbiome_16s_class_safe_colab.ipynb` as an Atacama-only 16S ASV notebook.
- Removed the synthetic `Soil_ASV_A` / `Soil_ASV_B` teaching path from the student-facing notebook.
- Removed visible form controls, live BLAST, neighbor joining, IQ-TREE, and bootstrap content from the student-facing notebook.
- Added `build_atacama_soil_asv_phylogeny_colab.py` and `verify_atacama_soil_asv_phylogeny_colab.py`.
- Added `classify_atacama_asvs_with_silva_static.py` and a cached SILVA nearest-reference taxonomy table.

## Data source

The generated cache uses the local QIIME 2 Atacama artifacts in `tmp/atacama_qiime2_source`:

- `atacama-table.qza`
- `atacama-rep-seqs.qza`
- `sample_metadata.tsv`

The generated table has 61 samples and 401 ASVs.

The builder now validates source artifacts before reading them:

- `.qza` files must be real zip archives and contain the expected QIIME payload (`feature-table.biom` or `dna-sequences.fasta`).
- metadata must contain the required Atacama columns, including `average-soil-relative-humidity`.
- the manifest records the local file path, source status, and SHA-256 hash for each input.
- broken 404 placeholders and non-QIIME files are rejected; the builder stops instead of synthesizing counts.

## Taxonomy cache

The final taxonomy cache now uses the ideal QIIME 2 route rather than the interim nearest-reference fallback.

Run environment:

- WSL2 Ubuntu imported onto `E:\WSL\qiime2-jammy`
- QIIME 2 Amplicon `2024.10.1`
- SILVA 138 99% OTUs Naive Bayes classifier for scikit-learn `1.4.2`
- classifier SHA256 verified as `c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616`

Command used inside WSL:

```text
qiime feature-classifier classify-sklearn \
  --i-classifier downloads/silva-138-99-nb-classifier.qza \
  --i-reads inputs/atacama-rep-seqs.qza \
  --o-classification outputs/atacama-taxonomy.qza \
  --p-n-jobs 1
```

Outputs:

- `E:\qiime2_atacama_taxonomy\outputs\atacama-taxonomy.qza`
- `E:\qiime2_atacama_taxonomy\outputs\exported-taxonomy\taxonomy.tsv`
- local copied source: `tmp\atacama_qiime2_source\atacama-taxonomy.tsv`
- published cache copy: `soil_16s_class_cache\goal2_atacama_qiime_taxonomy.tsv`

The older nearest-reference SILVA cache remains in the repo as a documented fallback only. The builder now prefers the QIIME `atacama-taxonomy.tsv`.

## Scientific honesty notes

- Taxonomic labels now come from QIIME 2 `feature-classifier classify-sklearn` output. They are still closest taxonomic matches, not species proof.
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
