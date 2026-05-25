from __future__ import annotations

import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "soil_microbiome_16s_class_safe_colab.ipynb"
OUTPUT_DIR = ROOT / "soil_microbiome_16s_outputs"
REPORT = ROOT / "soil_16s_class_cache" / "notebook_execution_report.json"


def main() -> None:
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)

    nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    namespace: dict[str, object] = {"__name__": "__notebook_validation__"}

    import matplotlib

    matplotlib.use("Agg")

    executed = 0
    for index, cell in enumerate(nb["cells"]):
        if cell.get("cell_type") != "code":
            continue
        source = cell.get("source", "")
        try:
            exec(compile(source, f"{NOTEBOOK.name}:cell_{index}", "exec"), namespace)
            executed += 1
            plt = namespace.get("plt")
            if plt is not None:
                plt.close("all")
        except Exception as exc:  # pragma: no cover - this is a validation utility
            raise RuntimeError(f"Notebook cell {index} failed") from exc

    expected_outputs = [
        OUTPUT_DIR / "soil_16s_distance_matrix.csv",
        OUTPUT_DIR / "soil_16s_closest_reference_report.csv",
        OUTPUT_DIR / "soil_16s_upgma_tree.newick",
        OUTPUT_DIR / "soil_16s_neighbor_joining_tree.newick",
        OUTPUT_DIR / "soil_16s_metadata_used.csv",
    ]
    missing = [str(path) for path in expected_outputs if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(f"Missing expected notebook outputs: {missing}")

    closest = namespace.get("closest")
    if closest is None:
        raise AssertionError("Notebook did not create closest-reference table")
    closest_map = dict(zip(closest["query"], closest["closest_reference"]))
    expected_closest = {
        "Soil_ASV_A": "Bacillus_subtilis_168",
        "Soil_ASV_B": "Rhizobium_leguminosarum_IAM12609",
    }
    if closest_map != expected_closest:
        raise AssertionError(f"Unexpected closest-reference map: {closest_map}")

    cached_blast_hits = namespace.get("cached_blast_hits")
    if cached_blast_hits is None or len(cached_blast_hits) != 10:
        raise AssertionError("Notebook did not parse the expected 10 cached BLAST-like XML hits")
    for required_column in ["bit_score", "teaching_e_value"]:
        if required_column not in cached_blast_hits.columns:
            raise AssertionError(f"Notebook did not parse cached BLAST-like XML column: {required_column}")

    report = {
        "status": "passed",
        "notebook": NOTEBOOK.name,
        "code_cells_executed": executed,
        "expected_outputs": [str(path.relative_to(ROOT)) for path in expected_outputs],
        "closest": closest_map,
        "cached_blast_xml_hits": int(len(cached_blast_hits)),
    }
    REPORT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
