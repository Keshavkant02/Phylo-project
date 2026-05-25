# Soil 16S Tree-Reading Printout: Answer Key

## Answers

1. `Soil_ASV_A` is closest to `Bacillus_subtilis_168`.

2. `Soil_ASV_B` is closest to `Rhizobium_leguminosarum_IAM12609`.

3. `Bacillus_subtilis_168` is not the ancestor of `Soil_ASV_A`. Both are tips. They are descendants of a shared ancestor in the tree.

4. C. Relatedness is decided by the most recent common ancestor.

5. No. Rotating a branch around a node changes the drawing, not the evolutionary relationship.

6. Accept either:
   - `Soil_ASV_A` and `Bacillus_subtilis_168`
   - `Soil_ASV_B` and `Rhizobium_leguminosarum_IAM12609`

7. No. High 16S identity supports a closest-reference claim, but one short 16S region does not prove exact species identity.

8. B. In this notebook, branch length is sequence-distance-derived branch length.

9. B. Bootstrap support is more appropriate for branch support. Adjusted p-values are for abundance or metadata association tests, not branch length.

10. Accept careful claims such as:

`Soil_ASV_A is closest to Bacillus_subtilis_168 in this cached 16S marker comparison, but this does not prove exact species identity.`

`Soil_ASV_B is closest to Rhizobium_leguminosarum_IAM12609 in this cached 16S marker comparison, but this does not prove exact species identity.`

## Instructor Notes

This activity is original and uses the soil 16S cache built for the Colab. It is inspired by the tree-thinking skill targets in Baum, Smith, and Donovan's "The Tree-Thinking Challenge" but does not reproduce the quiz pages.

Emphasize these misconceptions:

- Do not read trees as ladders of progress.
- Do not decide relatedness by the left-right order of tip labels.
- Do not treat living tips as ancestors of other living tips.
- Do not call branch length a p-value.
- Do not convert a closest-reference 16S hit into exact species proof.

Suggested class order:

1. Students answer the printout before running the notebook.
2. Students run the Colab and compare the actual UPGMA/NJ plots with the simplified tree.
3. Students revise one answer after seeing the distance matrix.
4. Students write one cautious report sentence.
