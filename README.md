# Atacama Soil 16S ASV Phylogeny Colab

Student-facing Colab workflow for reading a phylogenetic tree of real Atacama soil ASVs and connecting that tree to abundance, alpha diversity, and BH-corrected association tests.

Main notebook:

```text
soil_microbiome_16s_class_safe_colab.ipynb
```

The notebook is designed to run top-to-bottom in a vanilla Colab runtime. It uses real Atacama 16S ASVs from QIIME 2 Atacama tutorial artifacts and does not require live BLAST, a QIIME 2 install, or local software.

## Teaching Scope

- 16S amplicon microbiome data, not shotgun metagenomics.
- ASVs are exact cleaned 16S sequence patterns, not proved species.
- One distance-based UPGMA tree only.
- BH-adjusted q-values are for abundance-versus-metadata tests, not tree support.
- No synthetic `Soil_ASV_A` / `Soil_ASV_B` workflow in the student notebook.

## Build And Verify

Run:

```powershell
python build_atacama_soil_asv_phylogeny_colab.py
python verify_atacama_soil_asv_phylogeny_colab.py
```

Latest verification checks:

- 36 notebook cells.
- 46 QC-passed samples from 61 raw feature-table samples.
- 37 ASVs retained after the present-in-at-least-3-samples prevalence filter.
- 20 ASVs for relative abundance plots.
- 12 ASVs for the UPGMA tree.
- 8 ASVs for the alignment heatmap.
- 74 association-test rows: 37 ASVs x 2 metadata variables.
- 9 saved figure outputs for visual QA.

See `ATACAMA_GOAL2_QA_REPORT.md` and `soil_16s_class_cache/README.md` for provenance, cache dimensions, and scientific caveats.
