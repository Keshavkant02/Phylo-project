# Soil 16S Class Cache

This folder is the class-safe data cache for `soil_microbiome_16s_class_safe_colab.ipynb`.

The default classroom workflow should load these files instead of calling live BLAST, Entrez, SILVA, or other web services during class.

This cache supports a 16S marker/metabarcoding teaching workflow, not a shotgun metagenomics workflow.

One-command Colab pattern after this folder is pushed to GitHub:

```python
CACHE_BASE_URL = "https://raw.githubusercontent.com/<org>/<repo>/main/soil_16s_class_cache"
```

Files:

- `pilot_16s_references.fasta`: five NCBI 16S reference sequences.
- `pilot_16s_query_reads.fasta`: two synthetic teaching ASV/read sequences derived from the cached references.
- `pilot_16s_metadata.csv`: source, accession, phylum, soil context, and interpretation notes.
- `pilot_16s_cached_hits.csv`: precomputed closest-reference table for the teaching queries.
- `pilot_16s_cached_blast.xml`: cached BLAST-like XML for teaching ranked-hit parsing without live BLAST.
- `pilot_16s_abundance_table.csv`: toy microbiome-style counts.
- `pilot_16s_manifest.json`: provenance and class-use notes.
- `cache_validation_report.json`: last local validation result from the builder.
- `notebook_execution_report.json`: result from executing the generated notebook with the verifier.
- `VISUAL_QA_TUFTE.md`: visualization checklist used for the student-facing figures.
- `COLAB_ONE_CELL_LOADER.py`: copy-paste loader for any other Colab notebook.
- `PUBLISH_TO_GITHUB.md`: publish and pre-class verification checklist.

Retrieved date recorded in the manifest: 2026-05-25.
