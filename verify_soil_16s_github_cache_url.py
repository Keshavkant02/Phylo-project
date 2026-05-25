from __future__ import annotations

import argparse
import json
import urllib.request


CACHE_FILES = [
    "pilot_16s_references.fasta",
    "pilot_16s_query_reads.fasta",
    "pilot_16s_metadata.csv",
    "pilot_16s_cached_hits.csv",
    "pilot_16s_cached_blast.xml",
    "pilot_16s_abundance_table.csv",
    "pilot_16s_manifest.json",
    "cache_validation_report.json",
    "notebook_execution_report.json",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Verify that the soil 16S class cache is reachable from a raw GitHub base URL."
    )
    parser.add_argument(
        "base_url",
        help="Example: https://raw.githubusercontent.com/<org>/<repo>/main/soil_16s_class_cache",
    )
    args = parser.parse_args()

    base = args.base_url.rstrip("/")
    results = []
    for filename in CACHE_FILES:
        url = f"{base}/{filename}"
        with urllib.request.urlopen(url, timeout=20) as response:
            body = response.read()
        if not body:
            raise RuntimeError(f"{filename} downloaded as an empty file")
        results.append({"file": filename, "bytes": len(body)})

    print(json.dumps({"status": "passed", "base_url": base, "files": results}, indent=2))
    print()
    print("Use this in Colab:")
    print(f'CACHE_BASE_URL = "{base}"')


if __name__ == "__main__":
    main()
