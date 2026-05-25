# Atacama Soil Mini-Cache

This folder includes a small derived cache from the QIIME 2 Atacama soil microbiome tutorial.

Use this for browser-only teaching of real soil ASV abundance patterns. The heavy QIIME 2 steps are not run in Colab; this cache stores the derived tables students need.

Files:

- `atacama_sample_metadata_mini.csv`: selected sample metadata for table samples.
- `atacama_feature_table_top12.csv`: counts for the 12 most abundant ASVs plus all other ASVs.
- `atacama_relative_abundance_top12.csv`: relative abundance percentages for the same ASVs.
- `atacama_feature_key.csv`: ASV label to original QIIME feature ID mapping.
- `atacama_top_asv_sequences.fasta`: representative sequences for top ASVs.
- `atacama_top_asv_stats.csv`: Spearman humidity and Mann-Whitney vegetation tests with Benjamini-Hochberg q-values.
- `atacama_alpha_diversity.csv`: total reads, observed ASVs, and Shannon entropy per sample.
- `atacama_alpha_diversity_stats.csv`: alpha-diversity correlation with humidity and BH q-values.
- `atacama_mini_manifest.json`: provenance and source URLs.

Teaching caveat: adjusted p-values here test abundance/metadata associations. They are not tree branch support values.
