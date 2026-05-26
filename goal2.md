# Goal 2: Atacama-Only Soil ASV Phylogeny Colab

## Priority Order When Specs Conflict

When in doubt, prioritize in this order:

1. Scientific honesty
2. Pedagogical clarity
3. Visual polish
4. Strict adherence to specs below

If a spec rule conflicts with a clearer student outcome, break the spec rule and document the reason in a code comment.

## Target

- 60-75 minutes of class time.
- 28-36 notebook cells total.
- Run All works top-to-bottom on a vanilla Colab runtime.
- No student-facing form controls.
- No optional branches.
- No live BLAST, no QIIME 2 install inside Colab, no local software.

## Core Scientific Framing

- This is 16S amplicon microbiome data, not shotgun metagenomics.
- ASV = Amplicon Sequence Variant: an exact cleaned 16S sequence pattern after denoising.
- ASV is not the same as species. Use careful language: "ASV," "closest taxonomic match," "sequence relatedness." Never "proved species."
- BH-adjusted q-values apply to abundance/metadata association tests only.
- The phylogenetic tree is distance-based, not p-value-based.
- Use one UPGMA tree only. No neighbor joining, no IQ-TREE, no bootstrap.

## Required Section Order

Section 1: Story hook - Atacama Desert  
Section 2: Tree-thinking intro (mammals only)  
Section 3: Atacama dataset story  
Section 4: What is an ASV?  
Section 5: Load cached Atacama data  
Section 6: Alignment of representative ASV sequences  
Section 7: Distance matrix  
Section 8: UPGMA tree (the payoff)  
Section 9: Relative abundance - what is actually in these samples?  
Section 10: Alpha diversity  
Section 11: BH-corrected association tests  
Section 12: Final student report

The section numbering above is canonical. All verification, all cross-references, and all "see section N" markdown text must match these numbers exactly.

## Section 1 - Story Hook

Open with the Atacama Desert as the narrative anchor.

- Driest non-polar desert on Earth.
- Parts of it have measured zero rainfall for decades.
- Yet microbial life persists in the soil.
- Guiding question: Who lives there, where are they abundant, and how are they related?
- Keep to one short markdown cell. No code yet.

## Section 2 - Tree-Thinking Intro

Base on Baum, Smith and Donovan (2005), "The Tree-Thinking Challenge," *Science* 310:979-980.

Use familiar mammals only: dog, wolf, fox, bear, cat, lion. No microbes in this section.

Use exactly this Newick string for the mammal tree:

```text
(((dog:1,wolf:1):1,fox:2):1,(bear:3,(cat:2,lion:2):1.5):0.5);
```

Topology: dog and wolf are sister taxa; fox is sister to the dog/wolf clade; cat and lion are sister; bear is sister to the cat/lion clade.

Teach through short exercises that surface and correct misconceptions, not through bullet lists. Cover:

- Reading across tips: show same tree in three layouts (rectangular, ladderized, alternative tip order); ask which two are most related in each; reveal the answer is the same in all three.
- Node rotation: show same tree before and after rotating two internal nodes; ask "did relationships change?"; reveal: no.
- Tips vs inferred ancestors: one labeled diagram; tips = present-day organisms, internal nodes = inferred common ancestors whose DNA we do not have.
- Branch length axis: in this notebook, branch length = amount of DNA difference, not time. Make this explicit.
- Sister groups and MRCA: define both with a labeled diagram, demonstrate by tracing three pairs of taxa.

Budget: 4-6 cells total.

## Section 3 - Atacama Dataset Story

Explain that this is real QIIME 2 Atacama soil 16S data from the 10% subsample at:

```text
data.qiime2.org/2024.10/tutorials/atacama-soils/10p/
```

Briefly describe what was sampled: soil across a humidity/aridity gradient with vegetation metadata.

Keep technical provenance compact (one short paragraph). Cache mechanics belong in code comments or repo README, not in main teaching markdown.

## Section 4 - What Is An ASV?

High-school definition:

> An ASV is a precise DNA sequence pattern found after cleaning 16S sequencing reads.

Distinguish three different ideas that students will conflate if not addressed:

- Abundance: how much of an ASV is in a sample.
- Taxonomic match: what known group the ASV's sequence resembles.
- Evolutionary relatedness: how ASV sequences cluster in a tree.

These three concepts return in the synthesis question at the end. Plant the seed here.

## Section 5 - Load Cached Atacama Data

No visible controls. Fixed settings only. Load cached files silently from GitHub or embedded fallback.

Code comment must document the subset cascade:

```python
# ASV cascade - chosen for readability at each step:
#   top 50 by prevalence -> q-value tests (statistical power)
#   top 20 by mean abundance -> abundance plots (readable bars)
#   top 12 by mean abundance -> UPGMA tree (readable tip labels)
#   top 8  by mean abundance -> alignment heatmap (readable bases)
```

Cached parameters (do not deviate):

- Samples: 61
- ASVs in q-value tests: top 50 by prevalence (present in >= 10% of samples)
- ASVs in abundance plots: top 20 by mean relative abundance; all others collapsed to "Other"
- ASVs in UPGMA tree: top 12 by mean relative abundance
- ASVs in alignment heatmap: top 8 by mean relative abundance

Taxonomy assignments come from the SILVA-classified `taxonomy.tsv` distributed with the QIIME 2 Atacama tutorial. Display at genus level (or family if genus is unassigned). If neither is assigned, show "Unassigned at genus level" - do not invent taxonomy.

Student-facing output for this section is a small, polished dataset summary (sample count, ASV count, metadata columns shown), not a raw metadata dump.

## Section 6 - Alignment Of Representative ASV Sequences

Use the top 8 ASVs by mean relative abundance.

Explain alignment as: "putting DNA letters into comparable columns."

Alignment heatmap specification:

- 8 rows.
- 60-column window chosen for mixed conserved/variable content (not a gap-heavy region).
- Conserved columns (>= 87.5% identity across rows) rendered at 30% opacity.
- Variable columns rendered at 100% opacity.
- Compute per-column identity fraction, then multiply alpha. This is the single most important visual change - variable columns must visually carry the signal.
- Base colors (Okabe-Ito): A = `#56B4E9`, C = `#009E73`, G = `#E69F00`, T = `#CC79A7`.
- Gap `-` and N = `#EEEEEE` (very light gray).
- Remove: red mismatch outlines, top tick marks, row separator lines, bottom legend.
- Replace legend with a one-line caption below the figure: "Colored cells = DNA bases (A blue, C green, G orange, T pink). Faded columns = conserved across all sequences; bright columns = where the sequences differ."
- Figure: 12 inches wide; height = 0.5 inch per row + 1.2 inch margin.
- Y-axis: 10pt monospace, no underscores in displayed labels, right-aligned, about 1 cm breathing room before the heatmap.
- X-axis: tick only at first, middle, last column number.
- Title: "Aligned 16S marker window - variable columns carry the phylogenetic signal."
- Interpretation line below caption: "ASVs that look similar across the bright columns are more closely related. The phylogenetic tree (section 8) formalizes this intuition."

## Section 7 - Distance Matrix

Explain: "Small distance = more similar DNA sequence."

Build pairwise distances from the aligned ASV sequences (top 12 ASVs - same set used for the tree).

Plot a clean distance heatmap using viridis. No clutter, no oversized tick labels, no decorative styling. Annotate cell values if they fit cleanly; otherwise omit annotations.

## Section 8 - UPGMA Tree

Build one UPGMA tree only from the top 12 Atacama ASV representative sequences. No neighbor joining, no IQ-TREE, no bootstrap.

Plot a clean phylogram:

- Readable tip labels (ASV ID + closest genus in parentheses where available).
- Italicize genus names.
- Branch length axis labeled "branch length (sequence-distance units)".
- No top/right spines.
- No decorative boxes.
- No overlapping labels - if labels collide, increase figure height before reducing font size.

Explain in markdown:

> Tips are ASVs. Branch length shows sequence difference. Branches that meet recently (short path between them) suggest closer sequence relatedness - which is our best evidence of evolutionary relatedness, but not proof of it.

## Section 9 - Relative Abundance

Pivot the question: "What is actually in these soil samples?"

Use top 20 ASVs by mean relative abundance. Collapse all others into "Other."

Plot:

- Clean stacked bar or grouped bar showing relative abundance per sample (or by metadata category if more readable).
- Okabe-Ito categorical colors, cycled across the 20 ASVs (use a lookup so the same ASV gets the same color in the lollipop plot later).
- Direct labeling where possible; if a legend is needed, keep it compact and outside the plot area.

Table (styled with pandas.Styler, not raw dataframe):

- Columns: "ASV", "Closest taxonomic match", "Mean relative abundance (%)", "Maximum relative abundance (%)", "Samples detected (of 61)".
- Round all percentages to 2 decimals.
- Alternating row shading (`#fafafa` and white).

## Section 10 - Alpha Diversity

Explain: "Alpha diversity asks how diverse one sample is - how many ASV types are present, and how evenly distributed they are."

Show two metrics:

- Observed ASVs (number of ASV types per sample).
- Shannon diversity (richness plus evenness).

Plot diversity vs humidity (two side-by-side scatter plots, one per metric), with a faint trend line. Use Okabe-Ito blue for points, viridis for any continuous color overlay.

One interpretation line below the figure.

## Section 11 - BH-Corrected Association Tests

Pedagogy first. Open with this exact framing in a markdown cell:

> We tested 50 ASVs to see if abundance changes with humidity. With a p < 0.05 cutoff, we would expect about 50 x 0.05 = 2.5 false positives by random chance alone, even if nothing is truly associated. Benjamini-Hochberg (BH) correction adjusts for this. A q-value of 0.05 means we expect about 5% of discoveries below that threshold to be false alarms.

Statistical tests:

- Humidity (continuous): Spearman correlation between CLR-transformed ASV abundance and average-soil-relative-humidity.
- Vegetation (categorical, "yes" / "no"): Mann-Whitney U test on CLR-transformed abundance, two-sided.
- Multiple testing correction: BH FDR via `statsmodels.stats.multitest.multipletests(method="fdr_bh")`, applied separately per metadata variable across all 50 ASVs.

CLR transform implementation:

- Add a pseudo-count of 0.5 to all zero values before taking the log (standard microbiome convention to avoid log(0)).
- CLR per sample: `clr(x) = log(x) - mean(log(x))`.
- Document the pseudo-count choice in a code comment so students can see the assumption.

Computation note:

The q-value table is computed live in the notebook from the loaded ASV abundance table and sample metadata - not loaded from a precomputed CSV. "Cache" in this context means: compute once at notebook runtime, then reuse the resulting table across the q-value section's plots and tables. The computation must use real Atacama abundance values, never synthesized ones.

Plot 1: raw p vs q illustration

- 4 inch wide by 3 inch tall figure.
- Two side-by-side strip plots: all 50 raw p-values (left), all 50 q-values (right).
- Jittered points, Okabe-Ito blue, semi-transparent.
- Horizontal dashed line at y = 0.05 on both panels.
- Caption: "Multiple-testing correction is conservative on purpose - it is the cost of asking many questions at once."

Table: top 10 ASVs by smallest q-value, per metadata variable (one combined table with a metadata-variable column is acceptable to save cells).

Columns:

- "ASV"
- "Closest taxonomic match"
- "Mean abundance (%)"
- "Spearman rho" or "Effect size" (with arrow indicating direction)
- "Raw p-value" (scientific notation if < 0.001, else 3 decimals)
- "BH q-value" (same format)
- "Significant?" (green check if q < 0.05, light gray dash otherwise)

Style with pandas.Styler:

- Alternating row shading.
- Rounded values.
- q-value column with subtle gradient (light to slightly darker as q approaches 0).
- Significant rows marked with a thin `#009E73` left border.
- Width <= 1100px (no horizontal scrolling at 1100px viewport).

Plot 2: lollipop, not volcano

- Horizontal lollipop, one row per significant ASV (q < 0.05), sorted by effect size.
- X-axis: Spearman rho (humidity) or log2 fold-change (vegetation). One panel per metadata variable.
- X = 0 reference line.
- Dot size proportional to mean abundance.
- Dot color by phylum if available; otherwise Okabe-Ito blue.
- Y-axis labels: ASV ID + italicized genus in parentheses.
- Right-edge annotation: q-value in light gray.
- Caption: "Each bar is one ASV. Bars to the right are more abundant in higher-humidity samples; bars to the left, lower-humidity. Dot size = overall abundance. Only ASVs with BH q < 0.05 shown."

Do not use: volcano plots, log-scaled axes on the lollipop, significance stars (`*`, `**`, `***`).

## Section 12 - Final Student Report

Provide a structured student template in a markdown cell with `[your answer]` placeholders. Students answer:

1. Which 3 ASVs are most abundant in your samples? (refer to section 9 table)
2. Which ASVs are significantly associated with humidity? (refer to section 11 lollipop)
3. Which ASVs are significantly associated with vegetation? (refer to section 11 lollipop)
4. What does alpha diversity suggest about humid versus arid samples? (refer to section 10)
5. In the UPGMA tree, which ASVs cluster closest together? (refer to section 8)
6. Are closely related ASVs (section 8) also similar in abundance (section 9) or in humidity association (section 11)?
7. Synthesis: An ASV can be (a) very abundant, (b) statistically associated with humidity, and (c) closely related to another ASV in the tree - and these are three different things. Explain in 2-3 sentences why these are three different ideas and why all three matter.

The synthesis question (7) is the one assessment that tests whether the student actually understood the notebook.

## Cell Anatomy - Template For Every Results Section

Each results section follows this anatomy:

1. Markdown cell: 2-3 sentences explaining what students are about to see and why.
2. Code cell: loads data, computes, renders. One figure OR one styled table per code cell. No `print()` debugging, no `display(df.head())` leftovers.
3. Markdown cell: italic one-line interpretation of the result, plus an optional "Think:" prompt where relevant.

## Visual Design Rules

Before any figure, view the tufte skill (if available in the environment) and apply its matplotlib rcParams.

Every figure must have:

- Sentence-case title, left-aligned, slightly larger than axis labels.
- One-line italic caption below the figure with interpretation.
- No top/right spines.
- No gridlines unless they carry data.
- Direct labels preferred over legends.
- No 3D, gradients, drop shadows, decorative cards, or patterned fills.

Color rules:

- Categorical: Okabe-Ito.
  - DNA bases: A = `#56B4E9`, C = `#009E73`, G = `#E69F00`, T = `#CC79A7`.
  - Sample/category colors: `#0072B2`, `#E69F00`, `#009E73`, `#CC79A7`, `#56B4E9`, `#D55E00`, `#F0E442`, `#000000`.
- Continuous: viridis.
- Diverging: cividis or RdBu.
- Never use: jet, rainbow, hsv.

Tables:

- Always styled (pandas.Styler), never raw dataframe dumps.
- Rename technical column names into student-readable labels (for example, "Mean relative abundance (%)" not `mean_rel_abund_pct`).
- No horizontal scrolling at 1100px viewport width.
- Round numeric values appropriately for display.

Decorative boxes are allowed only to mark something pedagogically important: a Think prompt, a careful scientific claim, or an interpretation tip.

## Common Failure Modes To Avoid

The previous version had these failures. Explicitly avoid:

1. Inventing references when real data exists. Always use real Atacama ASVs and SILVA-derived taxonomy.
2. Narrating the pipeline instead of telling the science story. Markdown should be biology-first, not procedure-first.
3. Skipping interpretation lines. Every figure has one italic interpretation sentence. Non-negotiable.
4. Leaving debug `print()` calls and `display(df.head())` in cells. Strip these before final output.
5. Reusing one color scheme for both categorical and continuous data. Categorical = Okabe-Ito, continuous = viridis.
6. Fabricating statistics. Compute q-values live from the loaded abundance table; never invent p-values.
7. Letting taxonomy be guessed by the language model. If SILVA does not assign a name, the table shows "Unassigned at genus level."

## Remove From Student-Facing Notebook

- Synthetic `Soil_ASV_A` / `Soil_ASV_B` workflow.
- Five-reference NCBI controlled teaching layer.
- "Why cached data?" as a main lesson.
- Class-safe, pilot, controlled-teaching language anywhere in markdown.
- Visible `CACHE_BASE_URL` or `USE_GITHUB_CACHE` controls.
- `QUERY_TO_REPORT`, `TREE_METHOD_TO_SHOW`, sliders, or any form controls.
- Neighbor joining, IQ-TREE, bootstrap.
- BLAST-like teaching E-values.
- Toy sample plots.
- Bland metadata tables.

## Keep In Repo/Docs Only, Not In Student Notebook

- Cache provenance and source URLs.
- Builder scripts.
- Verifier scripts.
- Scientific QA report.
- Raw cache validation.

## Verification

Regenerate notebook from the builder script. Run notebook top-to-bottom with the verifier (`nbclient` or equivalent).

Confirm:

- 28-36 code+markdown cells total.
- 0 code cell errors.
- 0 visible form controls.
- 0 synthetic teaching sections remain.

Confirm outputs include:

- Styled dataset summary (section 5).
- Alignment heatmap with visible conserved-column desaturation (section 6).
- Distance matrix heatmap (section 7).
- UPGMA tree phylogram (section 8).
- Relative abundance styled table + plot (section 9).
- Alpha diversity table + scatter plots (section 10).
- Raw p vs q illustration plot (section 11).
- q-value styled results table (section 11).
- Lollipop association plot (section 11).
- Final report markdown template with all 7 placeholders (section 12).

Screenshot every figure cell and verify:

- No label overlap.
- No table horizontal scrolling at 1100px viewport width.
- Palettes match spec (spot-check 3 figures).
- Conserved-column desaturation is visibly different from variable columns in the alignment heatmap.
- Lollipop plot is sorted by effect size and directly labeled.
- Final report template has all 7 placeholders.

Push updated notebook, cache, and docs to GitHub.

## Compact Workflow Version

```text
Atacama story
  -> tree-reading with mammals
  -> real Atacama dataset
  -> ASVs
  -> load data
  -> alignment
  -> distance matrix
  -> UPGMA tree
  -> relative abundance
  -> alpha diversity
  -> BH q-values
  -> final student report
```
