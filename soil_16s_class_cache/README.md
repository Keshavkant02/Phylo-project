# Soil 16S Atacama ASV Cache

This folder supports the Atacama-only student notebook:

```text
soil_microbiome_16s_class_safe_colab.ipynb
```

The current student-facing notebook uses real Atacama 16S ASVs only. It does not use the older synthetic
`Soil_ASV_A` / `Soil_ASV_B` teaching workflow, live BLAST, neighbor joining, IQ-TREE, or bootstrap sections.

## Goal 2 files

- `goal2_atacama_sample_metadata.csv`: 46 QC-passed Atacama soil samples with humidity and vegetation metadata.
- `goal2_atacama_counts_retained_asvs.csv`: counts for the 37 ASVs present in at least 3 QC-passed samples, used for q-value tests.
- `goal2_atacama_relative_abundance_top20.csv`: top 20 retained ASVs by mean relative abundance, with remaining ASVs collapsed to `Other`.
- `goal2_atacama_alpha_diversity.csv`: observed ASVs and Shannon diversity per QC-passed sample.
- `goal2_atacama_feature_key.csv`: retained ASV labels, source feature IDs, abundance ranks, prevalence, and taxonomy fields.
- `goal2_atacama_rep_seqs_retained_asvs.fasta`: representative 16S sequences used for alignment, distance, and tree sections.
- `goal2_atacama_manifest.json`: source and cache-generation metadata.
- `goal2_verification_report.json`: latest verifier output.
- `goal2_runtime_audit.json`: runtime dimensions from the executed notebook.
- `goal2_figure_checks/`: saved PNG outputs from figure cells for visual QA.
- `ATACAMA_GOAL2_README.md`: provenance and subset policy.

## Scientific notes

Taxonomic labels come from the QIIME 2 SILVA taxonomy artifact where available. Labels are closest taxonomic
matches, not species proof. If SILVA does not provide a useful genus label, the notebook displays the most
specific honest rank available, such as family, order, class, phylum, or `Unassigned`.

The raw source table has 401 ASVs across 61 samples. The teaching notebook keeps 46 samples after read-depth
and metadata QC, then keeps 37 ASVs present in at least three of those samples. The sparseness is intentional:
it is part of the lesson, not hidden.
