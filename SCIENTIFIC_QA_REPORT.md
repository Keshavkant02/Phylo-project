# Scientific QA Report: Soil 16S Colab Pilot

Date checked: 2026-05-25

## Verdict

The notebook is suitable as a class-safe browser workflow for an introductory soil 16S phylogeny lesson. Students do not need to install local software. The default path uses an embedded cache, can optionally use raw GitHub cache files, and does not depend on live NCBI BLAST, Entrez, SILVA, or RDP calls during class.

Scientific scope must stay clear: this is a 16S marker/metabarcoding teaching workflow, not shotgun metagenomics, not a species-proof workflow, and not a publication-grade phylogenetic inference pipeline.

## Execution Evidence

Local notebook execution passed after rebuilding the notebook:

- `python build_soil_microbiome_16s_class_safe_pilot.py`
- `python verify_soil_microbiome_16s_class_safe_pilot.py`
- `python -m py_compile build_soil_microbiome_16s_class_safe_pilot.py verify_soil_microbiome_16s_class_safe_pilot.py verify_soil_16s_github_cache_url.py`

Observed validation result:

- 22 code cells executed.
- 5 reference 16S sequences loaded.
- 2 teaching ASV/query reads loaded.
- 10 cached BLAST-like XML hit rows parsed.
- 61 Atacama soil sample metadata rows loaded.
- 12 Atacama ASV statistic rows loaded.
- 3 Atacama alpha-diversity statistic rows loaded.
- 5 bootstrap-supported clades summarized from resampled alignment columns.
- `Soil_ASV_A -> Bacillus_subtilis_168`.
- `Soil_ASV_B -> Rhizobium_leguminosarum_IAM12609`.
- Distance matrix, closest-reference report, UPGMA Newick tree, NJ Newick tree, bootstrap support table, metadata export, Atacama ASV stats, Atacama alpha-diversity stats, Atacama metadata, and Atacama relative-abundance exports were written.

Raw GitHub cache check passed for the already-published cache URL:

```python
CACHE_BASE_URL = "https://raw.githubusercontent.com/Keshavkant02/Phylo-project/main/soil_16s_class_cache"
```

## Scientific Review

The workflow is scientifically honest for teaching because it separates three ideas that students often confuse:

1. A microbiome count table says a read/ASV exists in samples.
2. A marker comparison can identify the closest cached reference sequence.
3. A distance tree is a hypothesis from one aligned marker window, not proof of exact species identity.
4. Abundance statistics test sample-metadata associations and must not be read as tree support.

The notebook now states that this is a browser-only 16S marker/metabarcoding pilot and not shotgun metagenomics. The final reporting sentence is appropriately cautious: students report "closest reference" rather than "this is species X."

UPGMA and neighbor joining are appropriate for the first lesson because they make the distance-matrix-to-tree step visible. UPGMA carries an equal-rate/ultrametric assumption, so the notebook compares it with neighbor joining. For stronger later work, the notebook correctly points students toward MAFFT for multiple alignment and IQ-TREE for model-based maximum-likelihood inference.

## Code Review

Change made in this QA pass:

- Replaced the older `Bio.pairwise2` alignment path with Biopython's modern `Bio.Align.PairwiseAligner`.
- Added explicit 16S-vs-shotgun scope language to the generated notebook.
- Added an A/C/G/T/N/gap legend to the alignment-window plot.
- Kept student-facing identity percentages tied to the cached direct-hit table while using the shared distance matrix for tree construction.
- Rebuilt the notebook and cache reports.

Remaining deliberate simplifications:

- The cached `pilot_16s_cached_blast.xml` is BLAST-like teaching XML, not a live NCBI BLAST result.
- The star alignment is a class-safe teaching approximation. It is useful for a small reference/query set but should be replaced by MAFFT or another multiple aligner for project-grade analysis.
- The toy abundance table is for narrative context only. Real abundance/statistics teaching now uses the cached Atacama mini-dataset.
- The Atacama q-values are Benjamini-Hochberg adjusted association-test p-values. They are not bootstrap values and not phylogenetic branch support.
- The bootstrap-support table is a simple teaching resampling summary for the small UPGMA tree, not a publication-grade phylogenetic support analysis.
- No 3D effects are used in scientific figures because they would distort distances and reduce graphical integrity.

## Visualization Review

The figures follow the intended Tufte-style constraints:

- Okabe-Ito colors for categorical marks.
- `viridis` for continuous distance values.
- No decorative gradients, 3D charts, glow, or heavy borders.
- Direct axes and colorbar labels.
- Small alignment view with variable-column ticks and a compact base legend.

Before teaching live, run the notebook once in Colab at normal classroom projector width and inspect for label collisions in the tree plots and heatmap tick labels.

## Literature And Resource Backing

- Biopython `PairwiseAligner` supports global and local pairwise alignments and exposes alignment parameters: https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html
- Biopython `DistanceTreeConstructor` supports distance-matrix UPGMA and neighbor-joining tree construction: https://biopython.org/docs/latest/api/Bio.Phylo.TreeConstruction.html
- QIIME 2 Atacama soil tutorial provides the soil microbiome teaching dataset and sample metadata: https://docs.qiime2.org/2024.10/tutorials/atacama-soils/
- QIIME 2 q2-vsearch chimera tutorial provides the cached Atacama feature table and representative sequence artifacts used for the mini-cache: https://docs.qiime2.org/2024.10/tutorials/chimera/
- NCBI BLAST Common URL API supports programmatic submission/retrieval, but this is deliberately not the default class path: https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html
- SILVA provides aligned and quality-controlled SSU/LSU rRNA resources and web tools: https://www.arb-silva.de/documentation/background/
- RDP provides aligned and annotated rRNA sequence data plus classifier/aligner tools: https://academic.oup.com/nar/article/42/D1/D633/1063201
- GTDB is genome-based and current as Release 11-RS232 from 2026-04-15, useful as a taxonomy reference for advanced extension: https://gtdb.ecogenomic.org/
- 16S can often support genus-level identification but has known species-resolution limits; Janda and Abbott's review remains a standard caution: https://journals.asm.org/doi/10.1128/jcm.01228-07
- MAFFT is the right next step for real multiple sequence alignment work: https://pmc.ncbi.nlm.nih.gov/articles/PMC3603318/
- IQ-TREE is an appropriate model-based maximum-likelihood extension, with a web server and current releases: https://iqtree.github.io/
- Neighbor joining originates with Saitou and Nei's distance-based method: https://pubmed.ncbi.nlm.nih.gov/3447015/

## Go/No-Go

Go for a browser-only class pilot.

Do not present the output as exact species discovery. Present it as a controlled, reproducible closest-reference and relatedness exercise that prepares students for a stronger MAFFT/IQ-TREE or database-backed project workflow.
