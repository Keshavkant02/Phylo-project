"""One-cell loader for the class-safe soil 16S cache.

Paste this cell near the top of any Google Colab notebook after the
`soil_16s_class_cache` folder has been pushed to GitHub.

Edit CACHE_BASE_URL to your raw GitHub URL, then run the cell. It downloads
the cache files into `/content/soil_16s_class_cache`.
"""

from pathlib import Path
import urllib.request


CACHE_BASE_URL = "https://raw.githubusercontent.com/<org>/<repo>/main/soil_16s_class_cache"
CACHE_DIR = Path("/content/soil_16s_class_cache")
CACHE_FILES = [
    "pilot_16s_references.fasta",
    "pilot_16s_query_reads.fasta",
    "pilot_16s_metadata.csv",
    "pilot_16s_cached_hits.csv",
    "pilot_16s_cached_blast.xml",
    "pilot_16s_abundance_table.csv",
    "pilot_16s_manifest.json",
    "atacama_sample_metadata_mini.csv",
    "atacama_feature_table_top12.csv",
    "atacama_relative_abundance_top12.csv",
    "atacama_feature_key.csv",
    "atacama_top_asv_sequences.fasta",
    "atacama_top_asv_stats.csv",
    "atacama_alpha_diversity.csv",
    "atacama_alpha_diversity_stats.csv",
    "atacama_mini_manifest.json",
    "ATACAMA_MINI_README.md",
    "cache_validation_report.json",
    "notebook_execution_report.json",
]


CACHE_DIR.mkdir(exist_ok=True)
for filename in CACHE_FILES:
    url = f"{CACHE_BASE_URL.rstrip('/')}/{filename}"
    target = CACHE_DIR / filename
    with urllib.request.urlopen(url, timeout=30) as response:
        target.write_bytes(response.read())
    print(f"loaded {filename} -> {target}")

print(f"cache ready: {CACHE_DIR}")
