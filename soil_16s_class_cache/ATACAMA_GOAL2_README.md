# Goal 2 Atacama Soil ASV Cache

This cache supports the Atacama-only student Colab in `soil_microbiome_16s_class_safe_colab.ipynb`.
It is derived from real QIIME 2 2024.10 Atacama tutorial artifacts and contains no synthetic ASVs.

Source files used locally:

- a validated QIIME 2 feature table artifact containing `feature-table.biom`
- a validated QIIME 2 representative-sequence artifact containing `dna-sequences.fasta`
- Atacama sample metadata containing `average-soil-relative-humidity`

The builder can use already-local files in `tmp/atacama_qiime2_source/`, or fetch the
official QIIME tutorial artifacts when the network is available. Broken 404 files and
non-zip `.qza` placeholders are rejected. If real artifacts cannot be found, the builder
stops instead of synthesizing counts.

The notebook applies sample quality control first:

- raw output: 401 ASVs across 61 samples
- keep samples with at least 100 reads and complete humidity/vegetation metadata: 46 samples
- keep ASVs present in at least three QC-passed samples: 37 ASVs

The notebook then uses fixed subsets for readability:

- all 37 retained ASVs for association tests
- top 20 retained ASVs by mean relative abundance for abundance plots
- top 12 retained ASVs by mean relative abundance for the UPGMA tree
- top 8 retained ASVs by mean relative abundance for the alignment heatmap

Sparseness policy: the 10% tutorial subsample is shallow. The notebook does not hide that;
it makes quality control and conservative statistical claims part of the lesson.

Taxonomy policy: the preferred student-facing cache reads
`goal2_atacama_qiime_taxonomy.tsv`, produced by QIIME 2 `feature-classifier classify-sklearn`
with the SILVA 138 Naive Bayes classifier. The builder can also read local QIIME taxonomy
artifacts from `tmp/atacama_qiime2_source/`. Taxonomy remains a closest-match label, not species proof.
