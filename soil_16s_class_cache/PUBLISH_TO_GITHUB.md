# Publish The Soil 16S Class Cache

This folder is ready to be committed to any GitHub repository and loaded by Colab through raw GitHub URLs.

## Files To Commit

- `pilot_16s_references.fasta`
- `pilot_16s_query_reads.fasta`
- `pilot_16s_metadata.csv`
- `pilot_16s_cached_hits.csv`
- `pilot_16s_cached_blast.xml`
- `pilot_16s_abundance_table.csv`
- `pilot_16s_manifest.json`
- `cache_validation_report.json`
- `notebook_execution_report.json`
- `VISUAL_QA_TUFTE.md`
- `README.md`
- `COLAB_ONE_CELL_LOADER.py`

## Raw URL Pattern

After pushing this folder, set:

```python
CACHE_BASE_URL = "https://raw.githubusercontent.com/<org>/<repo>/main/soil_16s_class_cache"
```

If your default branch is `master`, use:

```python
CACHE_BASE_URL = "https://raw.githubusercontent.com/<org>/<repo>/master/soil_16s_class_cache"
```

## One-Cell Colab Loader

Paste the contents of `COLAB_ONE_CELL_LOADER.py` into any Colab notebook and replace `<org>/<repo>`.

The generated class notebook already has a safer version of this logic:

- `USE_GITHUB_CACHE=False` by default.
- If `USE_GITHUB_CACHE=True`, it tries the GitHub raw cache.
- If GitHub fails, it falls back to the embedded copy in the notebook.

## Pre-Class Check

Before teaching, run:

```python
import urllib.request

base = "https://raw.githubusercontent.com/<org>/<repo>/main/soil_16s_class_cache"
for name in [
    "pilot_16s_references.fasta",
    "pilot_16s_query_reads.fasta",
    "pilot_16s_metadata.csv",
    "pilot_16s_cached_hits.csv",
    "pilot_16s_cached_blast.xml",
    "pilot_16s_abundance_table.csv",
    "pilot_16s_manifest.json",
]:
    url = f"{base}/{name}"
    with urllib.request.urlopen(url, timeout=20) as r:
        print(name, len(r.read()), "bytes")
```

If that succeeds, the GitHub cache is reachable. If it fails during class, the generated notebook still has an embedded fallback.
