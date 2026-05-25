# Soil 16S Teaching Synthesis And Next Dataset Plan

Date: 2026-05-25

## Short Answer

The current Colab cache was built as a controlled teaching cache, not as a full real metagenomics dataset. That was a good decision for the first student run because it gives known answers, avoids live BLAST failures, and lets students learn the logic of species-like identification and tree reading without software installation.

For the fuller module you want, the design is now:

1. Keep the current tiny cache as the first "how the method works" lesson.
2. Use the cached QIIME 2 Atacama soil microbiome mini-dataset for real abundance/statistics.
3. Add a Baum-style tree-thinking handout using original questions based on our soil bacteria, not copied quiz pages.
4. Add a tree/evolution section that distinguishes closest-reference identification from evolutionary relatedness.

## How The Current Cache Was Built

The builder script is `build_soil_microbiome_16s_class_safe_pilot.py`.

The cache is generated in five steps:

1. Five known bacterial 16S rRNA sequences were embedded as FASTA records.
2. Each sequence was labeled with metadata: accession, organism name, phylum, source URL, soil context, and teaching note.
3. A shared 16S marker window was extracted using the common 16S start anchor `AGAGTTTGATCCTGGCTCAG`.
4. Two teaching ASV/query reads were created by taking two known marker windows and introducing a few substitutions.
5. Pairwise distances, cached ranked hits, and a BLAST-like XML file were written to `soil_16s_class_cache/`.

The important point to teach:

> We know the answer because we created the teaching queries from known references. That makes the first lesson robust and lets students focus on reasoning rather than debugging.

## Exact Cached Reference Set

All five reference records are from NCBI Nucleotide and were cached on 2026-05-25.

| Label | Organism | Accession | Why We Chose It |
|---|---|---|---|
| `Bacillus_subtilis_168` | *Bacillus subtilis* subsp. *subtilis* strain 168 | `NR_102783.2` | Common soil/rhizosphere model bacterium; Gram-positive decomposer-style example. |
| `Pseudomonas_fluorescens_CCM2115` | *Pseudomonas fluorescens* strain CCM 2115 | `NR_115715.1` | Rhizosphere-associated bacterium; useful plant-root contrast. |
| `Streptomyces_coelicolor_rrnD` | *Streptomyces coelicolor* | `Y00411.1` | Classic filamentous soil actinomycete; shows soil microbes are evolutionarily diverse. |
| `Rhizobium_leguminosarum_IAM12609` | *Rhizobium leguminosarum* type strain IAM 12609 | `D14513.1` | Root nodule/nitrogen-cycle example; includes ambiguous `N` bases from the source record. |
| `Acidobacterium_capsulatum_ATCC51196` | *Acidobacterium capsulatum* ATCC 51196 | `NR_074106.1` | Represents Acidobacteriota, a major soil-associated phylum. |

## Exact Teaching Queries

| Query | How It Was Made | Expected Closest Reference |
|---|---|---|
| `Soil_ASV_A` | Derived from the *Bacillus subtilis* 16S marker window with substitutions at positions 43, 141, 309, and 471. | `Bacillus_subtilis_168` |
| `Soil_ASV_B` | Derived from the *Rhizobium leguminosarum* 16S marker window with substitutions at positions 58, 220, 391, and 505. | `Rhizobium_leguminosarum_IAM12609` |

This means they are synthetic classroom reads, not newly discovered organisms.

## Do We Have E-Values?

Yes, but with an important caveat.

The cached BLAST-like XML includes `<Hsp_evalue>` and `<Hsp_bit-score>` fields. The builder assigns teaching values:

- Rank 1 hit: `1e-120`
- Rank 2 hit: `1e-80`
- Rank 3 hit: `1e-75`
- Rank 4 hit: `1e-70`
- Rank 5 hit: `1e-65`

These are not real NCBI BLAST E-values. They are teaching placeholders in a BLAST-like XML structure so students can learn what E-values and ranked hits look like without calling live BLAST.

Implemented in the current local notebook generator:

- Parse and show `E-value` and `bit score` in the hit table.
- Label them clearly as "cached teaching E-value" unless we later run real BLAST before class and cache its true XML output.

## Identification Versus Evolution

This is the key conceptual split for teaching.

### Species-Like Identification

Question:

> Which known reference sequence is my ASV closest to?

Evidence:

- Percent identity
- E-value or bit score if using BLAST-style search
- Ranked reference hits

Correct claim:

> This ASV is closest to a *Bacillus subtilis*-like cached reference.

Incorrect claim:

> This ASV proves we discovered *Bacillus subtilis* in the sample.

### Evolutionary Relatedness

Question:

> How are the references and queries related as a gene tree?

Evidence:

- Aligned 16S marker columns
- Pairwise distance matrix
- UPGMA or neighbor-joining tree
- Branch lengths in a phylogram
- Optional bootstrap support from resampled alignment columns

Correct claim:

> These sequences cluster by shared marker similarity. This tree is a 16S gene-tree hypothesis.

Incorrect claim:

> This one short 16S window proves the complete species tree.

## About P-Values, Adjusted P-Values, And Trees

Use different statistics for different questions.

For abundance questions:

- Example: does an ASV or phylum change across soil humidity, vegetation, depth, or site type?
- Use Mann-Whitney, Kruskal-Wallis, Spearman correlation, PERMANOVA, or similar tests depending on the question.
- If many ASVs/phyla are tested, use Benjamini-Hochberg false discovery rate adjustment.
- Report adjusted p-values or q-values.

For tree questions:

- Do not call branch lengths "p-values."
- A distance-based phylogram shows sequence distance as branch length.
- Use bootstrap support, not p-values, to show support for branches.
- For a simple Colab version, bootstrap by resampling alignment columns 100-500 times, rebuilding a tree, and counting how often a clade appears.

Recommended student language:

> Abundance statistics ask whether sample groups differ. Tree support asks whether a relationship is stable under resampling. They are related parts of the story, but not the same statistic.

## Real Soil Dataset Added

The current tiny teaching cache is good for a first controlled lesson. It is too small and synthetic for real abundance statistics, so the notebook now includes a second cached dataset for that purpose.

The added dataset is the QIIME 2 Atacama soil microbiome tutorial dataset.

Why it fits:

- It is real soil microbiome data.
- It is already used for teaching.
- It has environmental metadata: transect, pit, depth, elevation, and relative humidity.
- The biological story is clear: aridity/humidity is associated with soil microbiome structure.
- QIIME 2 already frames guiding questions for richness, evenness, beta diversity, metadata association, and phylum abundance.

How it is used in Colab:

- Do not run full QIIME 2 in the student Colab.
- Preprocess/cache small derived tables before class.
- Include:
  - sample metadata,
  - feature/ASV abundance table,
  - representative ASV sequences,
  - precomputed diversity/statistical summaries.

Implemented local cache files:

- `atacama_sample_metadata_mini.csv`: 61 Atacama soil samples with selected environmental metadata.
- `atacama_feature_table_top12.csv`: counts for the 12 most abundant ASVs plus all other ASVs.
- `atacama_relative_abundance_top12.csv`: relative abundance percentages.
- `atacama_feature_key.csv`: local ASV labels mapped to original QIIME feature IDs.
- `atacama_top_asv_sequences.fasta`: representative sequences for those top ASVs.
- `atacama_top_asv_stats.csv`: Spearman humidity tests and Mann-Whitney vegetation tests with Benjamini-Hochberg q-values.
- `atacama_alpha_diversity.csv`: total reads, observed ASVs, and Shannon entropy.
- `atacama_alpha_diversity_stats.csv`: alpha-diversity correlation with humidity and q-values.
- `atacama_mini_manifest.json`: source URLs, DOI, and provenance.

Students can now do real interpretation in the browser without installing conda, QIIME 2, or external binaries.

## Other Dataset Options

| Dataset | Good For | Concern |
|---|---|---|
| Current tiny cache | First tree/identity lesson | Synthetic, too small for real abundance statistics |
| QIIME 2 Moving Pictures | Canonical beginner microbiome tutorial | Human body-site data, not soil |
| QIIME 2 Atacama soil | Best fit for soil teaching | Need to cache derived tables so Colab stays simple |
| Earth Microbiome Project | Big, authentic global resource | Too large and complex for first class |
| USDA/plant rhizosphere tutorials | Nice plant-soil angle | Less clean as a stable browser-only class cache |

## Baum PDF: How To Use It Without Copying It

The local file `C:\Users\DELL\Downloads\baum.som.pdf` is the supporting online material for:

Baum, Smith, and Donovan, "The Tree-Thinking Challenge," *Science* 310, 979-980, 2005. DOI: `10.1126/science.1117727`.

The PDF contains tree-thinking quizzes and answer explanations. We should not copy the full worksheet into our materials. Instead, we should make our own printout using the same teaching ideas:

1. Closest relatives are determined by most recent common ancestor, not by which tips look physically close on the page.
2. Rotating branches around a node does not change the relationship.
3. Living tips are not ancestors of other living tips.
4. A clade contains an ancestor and all descendants.
5. Trait changes can be mapped onto a tree.
6. Bootstrap/support values are about confidence in branches, not abundance.
7. Different genes/windows can produce different trees.
8. A dendrogram of similarity is not automatically an evolutionary tree.

Implemented local handout artifacts:

- `soil_16s_tree_reading_printout.md`
- `soil_16s_tree_reading_answer_key.md`
- `soil_16s_tree_reading_printout.pdf`

These are original soil-bacteria questions inspired by Baum-style tree-thinking targets, not copied quiz pages.

## Proposed Student Module Flow

### Part 1: Read A Tree Before Running Code

Students learn:

- tip,
- node,
- branch,
- root,
- MRCA,
- clade,
- sister taxa,
- branch length,
- bootstrap/support.

Activity:

- Baum-style original soil-bacteria tree questions.

### Part 2: What Is In A Soil Sample?

Students load:

- toy or real cached abundance table,
- sample metadata.

Activity:

- Which ASVs are abundant in which samples?
- Are counts the same as species identity? No.

### Part 3: What Is This ASV Closest To?

Students use:

- cached reference sequences,
- cached ranked hit table,
- percent identity,
- E-value/bit score if shown.

Activity:

- Write a careful closest-reference claim.

### Part 4: How Are These Organisms Related?

Students use:

- marker window,
- alignment view,
- distance matrix,
- UPGMA and NJ trees.

Activity:

- Identify sister relationships and clades.
- Explain why tree topology matters more than left-right order.

### Part 5: Real Soil Abundance And Statistics

Use Atacama-derived cached data already embedded in the generated notebook.

Activities:

- Plot phylum abundance across humidity/aridity.
- Run simple tests or correlations.
- Adjust p-values with Benjamini-Hochberg.
- Interpret q-values cautiously.

### Part 6: Final Report

Students write:

1. what the ASV is closest to,
2. how the references are evolutionarily related,
3. whether abundance patterns differ across soil metadata,
4. what the analysis cannot prove.

## Remaining Next Implementation

1. Add a bootstrap-support teaching plot for the small 16S tree.
2. Add one short "read this tree" cell that directly references the printout concepts before students interpret the generated tree.
3. Keep all heavy preprocessing outside the student Colab.
4. Later option: replace teaching E-values with real pre-run BLAST XML if exact BLAST statistics are important for the class.

## Sources Reviewed

- Current local cache files under `soil_16s_class_cache/`.
- Builder script `build_soil_microbiome_16s_class_safe_pilot.py`.
- Local Baum supporting-material PDF: `C:\Users\DELL\Downloads\baum.som.pdf`.
- QIIME 2 Moving Pictures tutorial: https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html
- QIIME 2 Atacama soil tutorial: https://docs.qiime2.org/2024.10/tutorials/atacama-soils/
- Earth Microbiome Project: https://earthmicrobiome.org/
- MIT 6.047/6.878 molecular evolution and phylogenetics lecture: https://web.mit.edu/6.047/book-2012/Lecture20_Phylogenetics/Lecture20_Phylogenetics_standalone.pdf
- Baum, Smith, Donovan tree-thinking material: https://arboretum.harvard.edu/wp-content/uploads/2020/07/phylogeny-quizes.pdf
