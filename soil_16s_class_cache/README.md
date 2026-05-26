# Soil 16S Atacama ASV Cache

This folder supports the Atacama-only student notebook:

```text
soil_microbiome_16s_class_safe_colab.ipynb
```

The current student-facing notebook uses real Atacama 16S ASVs only. It does not use the older synthetic
`Soil_ASV_A` / `Soil_ASV_B` teaching workflow, live BLAST, neighbor joining, IQ-TREE, or bootstrap sections.

## Goal 2 files

- `goal2_atacama_sample_metadata.csv`: 61 Atacama soil samples with humidity and vegetation metadata.
- `goal2_atacama_counts_top50.csv`: counts for the top 50 ASVs by prevalence, used for q-value tests.
- `goal2_atacama_relative_abundance_top20.csv`: top 20 ASVs by mean relative abundance, with remaining ASVs collapsed to `Other`.
- `goal2_atacama_alpha_diversity.csv`: observed ASVs and Shannon diversity per sample.
- `goal2_atacama_feature_key.csv`: ASV labels, source feature IDs, abundance ranks, prevalence, and taxonomy fields.
- `goal2_atacama_rep_seqs_top50_union.fasta`: representative 16S sequences used for alignment, distance, and tree sections.
- `goal2_atacama_manifest.json`: source and cache-generation metadata.
- `goal2_verification_report.json`: latest verifier output.
- `goal2_runtime_audit.json`: runtime dimensions from the executed notebook.
- `goal2_figure_checks/`: saved PNG outputs from figure cells for visual QA.
- `ATACAMA_GOAL2_README.md`: provenance and subset policy.

## Scientific notes

Taxonomy is not inferred. If a valid SILVA taxonomy artifact is not present in the local source directory,
the notebook displays `Unassigned at genus level`.

The source table has only nine ASVs present in at least 10% of samples. The q-value section tests the requested
top 50 ASVs by prevalence, and the lollipop plot labels only discoveries that also meet the 10% prevalence
threshold.
