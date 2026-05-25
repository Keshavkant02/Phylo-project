# Soil 16S Phylogeny Class Pilot

Class-safe Colab workflow for a soil-microbiome-flavored introduction to 16S marker-gene phylogeny.

## Main Notebook

- `soil_microbiome_16s_class_safe_colab.ipynb`

Students can run the notebook top to bottom. The default mode uses an embedded cache, so live BLAST/NCBI calls are not required during class.

## GitHub Cache URL

After this repository is pushed to GitHub, use this in Colab:

```python
CACHE_BASE_URL = "https://raw.githubusercontent.com/Keshavkant02/Phylo-project/main/soil_16s_class_cache"
USE_GITHUB_CACHE = True
```

The notebook falls back to its embedded cache if the raw GitHub cache is unreachable.

## Cache Contents

- `soil_16s_class_cache/pilot_16s_references.fasta`
- `soil_16s_class_cache/pilot_16s_query_reads.fasta`
- `soil_16s_class_cache/pilot_16s_cached_hits.csv`
- `soil_16s_class_cache/pilot_16s_cached_blast.xml`
- `soil_16s_class_cache/pilot_16s_metadata.csv`
- `soil_16s_class_cache/pilot_16s_abundance_table.csv`

## Validation

Run:

```powershell
python verify_soil_microbiome_16s_class_safe_pilot.py
python verify_soil_16s_github_cache_url.py "https://raw.githubusercontent.com/Keshavkant02/Phylo-project/main/soil_16s_class_cache"
```

Expected closest-reference calls:

- `Soil_ASV_A -> Bacillus_subtilis_168`
- `Soil_ASV_B -> Rhizobium_leguminosarum_IAM12609`

