# Goal 2: Atacama Soil 16S ASV Teaching Notebook

Build and maintain a polished, single-path Google Colab notebook for high-school students using real QIIME 2 Atacama soil 16S data.

Keystone: do not hide the sparseness; make it part of the lesson.

## Required Workflow

1. Section 1: Story hook - Atacama Desert
2. Section 2: Tree-thinking intro (mammals only)
3. Section 3: Atacama dataset story (+ small provenance banner)
4. Section 4: What is an ASV?
5. Section 5: Load data and apply QC
6. Section 6: Alignment of representative ASV sequences
7. Section 7: Distance matrix
8. Section 8: UPGMA tree (the payoff)
9. Section 9: Relative abundance - what is actually in these samples?
10. Section 10: Alpha diversity
11. Section 11: BH-corrected association tests
12. Section 12: Final student report
13. Appendix: Full data provenance

## Scientific Rules

- This is 16S amplicon microbiome data, not shotgun metagenomics.
- ASV means Amplicon Sequence Variant: a precise cleaned 16S sequence pattern after DADA2 denoising.
- ASV is not the same as species. Use "ASV", "closest taxonomic match", and "sequence relatedness"; never "proved species".
- The tree is a distance-based UPGMA tree from sequence-distance units. Use no neighbor joining, IQ-TREE, or bootstrap.
- BH-adjusted q-values apply only to abundance/metadata association tests, not tree branch support.

## Real Data Dimensions

- Raw source: 401 ASVs across 61 samples.
- Sample QC: keep samples with at least 100 reads and complete humidity/vegetation metadata.
- Sample QC result: 46 samples.
- ASV prevalence filter: keep ASVs present in at least 3 of the 46 QC-passed samples.
- ASV filter result: 37 ASVs.

Cascade after filtering:

- q-value tests: all 37 retained ASVs.
- abundance plot/table: top 20 retained ASVs by mean relative abundance, with other ASVs collapsed to `Other`.
- UPGMA tree: top 12 retained ASVs by mean relative abundance.
- alignment heatmap: top 8 retained ASVs by mean relative abundance.

## Verification Target

Run:

```powershell
python build_atacama_soil_asv_phylogeny_colab.py
python verify_atacama_soil_asv_phylogeny_colab.py
```

Verifier must confirm:

- 36 notebook cells.
- 0 code-cell execution errors.
- 46 QC-passed samples.
- 37 retained ASVs.
- 74 association-test rows.
- final report has 7 `[your answer]` placeholders.
- no synthetic teaching workflow, form controls, BLAST-like E-values, neighbor joining, IQ-TREE, or bootstrap.
