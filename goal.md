# Goal: Class-Safe Soil 16S Phylogeny Pilot

## Objective

Build a polished Google Colab student app for a soil-microbiome-flavored introduction to phylogenetic tree analysis.

The class priority is reliability: students should open the Colab, run top to bottom, and complete the analysis without installing local software, waiting on live BLAST, or debugging web-service failures.

## Pilot Scope

- Theme: soil microbiome species finding and relatedness reporting.
- Dataset: five soil-relevant 16S rRNA reference species cached before class.
- Student task: compare one or more unknown soil ASV/read sequences against cached references, identify closest references, build distance-based trees, and report the claim cautiously.
- Teaching level: introductory but scientifically honest.
- Main method: cached 16S marker-window comparison, distance matrix, UPGMA, and neighbor joining.
- Optional extension: IQ-TREE for maximum-likelihood comparison after students understand distance trees.

## Class-Safe Data Strategy

Default class mode must never depend on live NCBI/BLAST calls.

Required cache files:

- `pilot_16s_references.fasta`: five vetted NCBI 16S reference sequences.
- `pilot_16s_query_reads.fasta`: one or more teaching ASV/read sequences derived from the reference set and clearly marked as teaching queries.
- `pilot_16s_metadata.csv`: labels, accessions, source database, source URL, phylum/group, soil context, and interpretation notes.
- `pilot_16s_abundance_table.csv`: small toy sample table so the workflow feels microbiome-like without pretending to be full shotgun metagenomics.
- `pilot_16s_manifest.json`: cache provenance, retrieval date, and class-use notes.

Colab loading rule:

```python
CACHE_BASE_URL = "https://raw.githubusercontent.com/<org>/<repo>/main/soil_16s_class_cache"
```

The notebook should try GitHub raw files when `USE_GITHUB_CACHE=True`, but it must also contain an embedded fallback copy of the same cache. That makes the class robust even if GitHub/raw access is slow during a session.

## Workflow Map

1. Start with the question: which known soil bacteria are our ASV/read sequences closest to?
2. Load cached 16S references and query reads.
3. Inspect metadata and source provenance.
4. Extract a shared 16S marker window for a controlled teaching comparison.
5. View a small alignment window and identify conserved/variable columns.
6. Compute pairwise sequence distances.
7. Visualize the distance matrix with `viridis` or `cividis`.
8. Build UPGMA and neighbor-joining trees with Biopython.
9. Report closest references and clades.
10. State what the tree can and cannot claim from one short 16S marker.
11. Optional: replace the cached teaching reads with project reads later.
12. Optional: run MAFFT/IQ-TREE as an advanced model-based extension.

## Tool Decision

Use Biopython for the main class path:

- `Bio.Align.PairwiseAligner` for the small browser-safe star alignment.
- `DistanceTreeConstructor.upgma` for UPGMA.
- `DistanceTreeConstructor.nj` for neighbor joining.
- No external binary is needed for the basic class-safe tree.

Use IQ-TREE only as an extension:

- IQ-TREE is excellent for model-based maximum-likelihood phylogenetics.
- It is not the simplest tool for teaching UPGMA.
- Once students understand distance trees, IQ-TREE can show how a model-based tree differs from UPGMA/NJ.

## Visualization Rules

Use Tufte-style discipline:

- Every visual mark must carry data or orientation.
- Use direct labels where possible.
- Avoid decorative cards, glow, heavy borders, and cramped legends.
- Use Okabe-Ito colors for categories.
- Use `viridis` or `cividis` for continuous distances.
- Keep labels outside branch geometry where possible.
- Apply the eraser test: remove anything that does not teach.
- Apply the collision test: no label should overlap a mark or another label.

Recommended palette:

- orange `#E69F00`
- sky blue `#56B4E9`
- bluish green `#009E73`
- yellow `#F0E442`
- blue `#0072B2`
- vermillion `#D55E00`
- reddish purple `#CC79A7`
- off-black `#222222`
- light guide gray `#DDDDDD`

## Student Reporting Target

Students should finish with a claim like:

> Soil_ASV_A clusters closest to the Bacillus subtilis-like reference in this cached 16S marker comparison. Because this is one short 16S region, I can report a closest reference, not prove exact species identity or a complete species tree.

## Acceptance Checks

- `goal.md` exists and records the class-safe pilot design.
- A new Colab notebook exists for the soil 16S pilot.
- Cache files exist under `soil_16s_class_cache/`.
- Notebook has Colab form controls.
- Notebook can load from GitHub raw cache or embedded fallback.
- Notebook avoids live BLAST/NCBI in the default path.
- Notebook includes UPGMA and neighbor-joining trees.
- Notebook explains that IQ-TREE is optional and model-based.
- Notebook code cells parse as Python.
- Plotting code uses Okabe-Ito and viridis/cividis.
- Claims are cautious about 16S species resolution.
- Scientific QA report records evidence, literature support, limitations, and go/no-go verdict.
