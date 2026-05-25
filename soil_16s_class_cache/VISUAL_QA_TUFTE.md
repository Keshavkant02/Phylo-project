# Visual QA: Tufte Checklist

This checklist applies to `soil_microbiome_16s_class_safe_colab.ipynb`.

## Palette

- Categories use Okabe-Ito colors.
- Continuous distance values use `viridis`.
- No decorative gradients, glow, or heavy card styling.

## Figure Checks

- Workflow map: labels sit above the marks; connectors are thin gray lines.
- Abundance plot: stacked bars encode the toy count table; no 3D effects.
- Atacama scatterplots: point positions encode real samples; vegetation is color-coded with Okabe-Ito colors.
- Atacama q-value plot: x-axis uses BH adjusted p-values and marks the 0.05 guide line.
- Alignment view: base colors encode A/C/G/T/N/gap and include a compact legend; variable columns are small ticks.
- Distance matrix: colorbar says exactly what the values mean.
- Tree plots: branch length axis remains visible; unnecessary plot borders are removed.

## Eraser Test

Every visible element should either encode data, orient the student, or label the data. Remove decorations that do none of these jobs.

## Collision Test

Before teaching, run the notebook in Colab at normal width and inspect:

- axis labels,
- tip labels,
- workflow labels,
- table wrapping,
- heatmap tick labels.

No label should overlap another label or hide a data mark.

## Scientific Integrity

- The distance heatmap uses a zero baseline and a monotonic sequential colormap.
- The notebook states that the tree is built from one 16S marker window.
- The report sentence says closest reference, not exact species proof.
