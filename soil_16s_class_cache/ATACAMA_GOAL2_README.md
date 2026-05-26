# Goal 2 Atacama Soil ASV Cache

This cache supports the Atacama-only student Colab in `soil_microbiome_16s_class_safe_colab.ipynb`.
It is derived from QIIME 2 2024.10 Atacama tutorial artifacts and contains no synthetic ASVs.

Source files used locally:

- `tmp/atacama_qiime2_source/atacama-table.qza`
- `tmp/atacama_qiime2_source/atacama-rep-seqs.qza`
- `tmp/atacama_qiime2_source/sample_metadata.tsv`

The notebook uses fixed subsets for readability:

- top 50 ASVs by prevalence for association tests
- top 20 ASVs by mean relative abundance for abundance plots
- top 12 ASVs by mean relative abundance for the UPGMA tree
- top 8 ASVs by mean relative abundance for the alignment heatmap

Prevalence policy: the source artifact has only nine ASVs present in at least 10% of samples,
so the notebook uses the top 50 by prevalence for the BH correction lesson and treats the
lowest-prevalence ASVs cautiously.

Taxonomy policy: if a real SILVA-style taxonomy artifact is available in the source directory, the
builder reads it. If not, taxonomy is shown as `Unassigned at genus level`; the builder does not
infer or invent taxonomic matches.
