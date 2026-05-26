# Goal 2 Atacama Soil ASV Cache

This cache supports the Atacama-only student Colab in `soil_microbiome_16s_class_safe_colab.ipynb`.
It is derived from QIIME 2 2024.10 Atacama tutorial artifacts and contains no synthetic ASVs.

Source files used locally:

- a validated QIIME 2 feature table artifact containing `feature-table.biom`
- a validated QIIME 2 representative-sequence artifact containing `dna-sequences.fasta`
- Atacama sample metadata containing `average-soil-relative-humidity`

The builder can use already-local files in `tmp/atacama_qiime2_source/`, or fetch the
official QIIME tutorial artifacts when the network is available. Broken 404 files and
non-zip `.qza` placeholders are rejected. If real artifacts cannot be found, the builder
stops instead of synthesizing counts.

The notebook uses fixed subsets for readability:

- top 50 ASVs by prevalence for association tests
- top 20 ASVs by mean relative abundance for abundance plots
- top 12 ASVs by mean relative abundance for the UPGMA tree
- top 8 ASVs by mean relative abundance for the alignment heatmap

Prevalence policy: the source artifact has only nine ASVs present in at least 10% of samples,
so the notebook uses the top 50 by prevalence for the BH correction lesson and treats the
lowest-prevalence ASVs cautiously.

Taxonomy policy: the preferred student-facing cache reads `goal2_atacama_silva_static_taxonomy_assignments.csv`,
which assigns each Atacama ASV to its nearest SILVA 138 515F/806R reference sequence by local alignment.
This is not QIIME 2 Naive Bayes classification, and it is not species proof. If that cache is absent,
the builder can read a real QIIME/SILVA taxonomy artifact from the source directory; otherwise taxonomy
is shown as `Unassigned at genus level`.
