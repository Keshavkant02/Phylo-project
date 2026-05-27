from __future__ import annotations

import base64
import json
import re
from pathlib import Path

import nbformat
from nbclient import NotebookClient


ROOT = Path(__file__).resolve().parent
NOTEBOOK_PATH = ROOT / "soil_microbiome_16s_class_safe_colab.ipynb"
QA_DIR = ROOT / "soil_16s_class_cache" / "goal2_figure_checks"
AUDIT_PATH = ROOT / "soil_16s_class_cache" / "goal2_verification_report.json"
EXECUTED_PATH = ROOT / "soil_16s_class_cache" / "goal2_executed_notebook.ipynb"


FORBIDDEN = [
    "Soil_ASV_A",
    "Soil_ASV_B",
    "USE_GITHUB_CACHE",
    "CACHE_BASE_URL",
    "QUERY_TO_REPORT",
    "TREE_METHOD_TO_SHOW",
    "#@param",
    "neighbor joining",
    "IQ-TREE",
    "bootstrap",
    "BLAST-like",
    "toy sample",
    "top50",
    "counts_top50",
    "in_q_value_top50",
]

EXPECTED_SECTIONS = [
    "Section 1: Story hook - Atacama Desert",
    "Section 2: Tree-thinking intro (mammals only)",
    "Section 3: Atacama dataset story",
    "Section 4: What is an ASV?",
    "Section 5: Load data and apply QC",
    "Section 6: Alignment of representative ASV sequences",
    "Section 7: Distance matrix",
    "Section 8: UPGMA tree (the payoff)",
    "Section 9: Relative abundance - what is actually in these samples?",
    "Section 10: Alpha diversity",
    "Section 11: BH-corrected association tests",
    "Section 12: Final student report",
]


def cell_source(cell) -> str:
    source = cell.get("source", "")
    if isinstance(source, list):
        return "".join(source)
    return source


def verify_static(nb) -> dict[str, object]:
    all_source = "\n".join(cell_source(cell) for cell in nb.cells)
    markdown_source = "\n".join(cell_source(cell) for cell in nb.cells if cell.get("cell_type") == "markdown")
    code_source = "\n".join(cell_source(cell) for cell in nb.cells if cell.get("cell_type") == "code")

    cell_count = len(nb.cells)
    if not (28 <= cell_count <= 36):
        raise AssertionError(f"Expected 28-36 cells, found {cell_count}")

    found_bad = [needle for needle in FORBIDDEN if needle.lower() in all_source.lower()]
    if found_bad:
        raise AssertionError(f"Forbidden old workflow text remains: {found_bad}")

    sections = re.findall(r"Section \d+: [^\n]+", markdown_source)
    if sections != EXPECTED_SECTIONS:
        raise AssertionError(f"Section order mismatch:\n{sections}")

    if "[your answer]" not in markdown_source or markdown_source.count("[your answer]") != 7:
        raise AssertionError("Final report must contain exactly seven [your answer] placeholders")

    if "37 × 0.05 ≈ 2" not in markdown_source:
        raise AssertionError("Multiple-testing explanation must use 37 × 0.05 ≈ 2")
    if "401 ASVs across 61 samples in the raw output; 37 ASVs across 46 samples after QC" not in markdown_source:
        raise AssertionError("Section 3 provenance banner is missing or has the wrong counts")
    if "Appendix: Full data provenance" not in markdown_source:
        raise AssertionError("Technical provenance appendix is missing")

    if "print(" in code_source:
        raise AssertionError("Student notebook contains print() calls")

    return {
        "cell_count": cell_count,
        "code_cells": sum(cell.get("cell_type") == "code" for cell in nb.cells),
        "markdown_cells": sum(cell.get("cell_type") == "markdown" for cell in nb.cells),
        "sections": sections,
    }


def append_runtime_audit(nb) -> None:
    audit_code = r"""
import json
from pathlib import Path

runtime_audit = {
    "samples": int(len(metadata)),
    "q_value_asvs": int(len(q_asvs)),
    "abundance_asvs": int(len(top20_asvs)),
    "tree_asvs": int(len(top12_asvs)),
    "alignment_asvs": int(len(top8_asvs)),
    "association_rows": int(len(association_results)),
    "association_variables": sorted(association_results["Variable"].unique().tolist()),
    "significant_rows": int((association_results["BH q-value"] < 0.05).sum()),
    "distance_matrix_shape": list(distance_matrix.shape),
    "tree_newick_endswith_semicolon": bool(tree_newick.endswith(";")),
    "taxonomy_note": manifest["taxonomy_note"],
    "manifest_retained_asvs": int(manifest["retained_asvs"]),
    "samples_dropped_by_qc": int(manifest["samples_dropped_by_qc"]),
    "min_total_reads": int(metadata["total_reads"].min()),
}
with Path("soil_16s_class_cache/goal2_runtime_audit.json").open("w", encoding="utf-8", newline="\n") as handle:
    handle.write(json.dumps(runtime_audit, indent=2) + "\n")
"""
    nb.cells.append(nbformat.v4.new_code_cell(audit_code))


def execute_notebook(nb):
    client = NotebookClient(nb, timeout=900, kernel_name="python3", resources={"metadata": {"path": str(ROOT)}})
    return client.execute()


def extract_outputs(nb) -> dict[str, object]:
    QA_DIR.mkdir(parents=True, exist_ok=True)
    for old_png in QA_DIR.glob("figure_cell*_output*.png"):
        old_png.unlink()

    figure_files = []
    html_tables = 0
    for cell_index, cell in enumerate(nb.cells):
        if cell.get("cell_type") != "code":
            continue
        for output_index, output in enumerate(cell.get("outputs", [])):
            data = output.get("data", {}) if isinstance(output, dict) else {}
            if "image/png" in data:
                raw = data["image/png"]
                if isinstance(raw, list):
                    raw = "".join(raw)
                out_path = QA_DIR / f"figure_cell{cell_index:02d}_output{output_index:02d}.png"
                out_path.write_bytes(base64.b64decode(raw))
                figure_files.append(str(out_path.relative_to(ROOT)))
            if "text/html" in data:
                html = data["text/html"]
                if isinstance(html, list):
                    html = "".join(html)
                if "<table" in html:
                    html_tables += 1
                    widths = [int(value) for value in re.findall(r"width:\s*(\d+)px", html)]
                    if widths and max(widths) > 1100:
                        raise AssertionError(f"Table width exceeds 1100px: {max(widths)}")
            text = output.get("text", "")
            if isinstance(text, list):
                text = "".join(text)
            if re.search(r"divide by zero|RuntimeWarning", text, flags=re.IGNORECASE):
                raise AssertionError(f"Runtime warning leaked into notebook output: {text[:200]}")

    if len(figure_files) < 8:
        raise AssertionError(f"Expected at least 8 figure screenshots, found {len(figure_files)}")
    if html_tables < 4:
        raise AssertionError(f"Expected at least 4 styled table outputs, found {html_tables}")

    return {"figure_files": figure_files, "styled_table_outputs": html_tables}


def verify_runtime_audit() -> dict[str, object]:
    path = ROOT / "soil_16s_class_cache" / "goal2_runtime_audit.json"
    audit = json.loads(path.read_text(encoding="utf-8"))
    expected = {
        "samples": 46,
        "q_value_asvs": 37,
        "abundance_asvs": 20,
        "tree_asvs": 12,
        "alignment_asvs": 8,
        "association_rows": 74,
        "distance_matrix_shape": [12, 12],
        "manifest_retained_asvs": 37,
        "samples_dropped_by_qc": 15,
        "min_total_reads": 138,
    }
    for key, value in expected.items():
        if audit.get(key) != value:
            raise AssertionError(f"Runtime audit mismatch for {key}: expected {value}, found {audit.get(key)}")
    if audit.get("association_variables") != ["Humidity", "Vegetation"]:
        raise AssertionError(f"Unexpected association variables: {audit.get('association_variables')}")
    if not audit.get("tree_newick_endswith_semicolon"):
        raise AssertionError("UPGMA Newick string was not created")
    return audit


def main() -> None:
    nb = nbformat.read(NOTEBOOK_PATH, as_version=4)
    static_report = verify_static(nb)
    append_runtime_audit(nb)
    executed = execute_notebook(nb)
    nbformat.write(executed, EXECUTED_PATH)
    output_report = extract_outputs(executed)
    runtime_report = verify_runtime_audit()

    report = {
        "notebook": str(NOTEBOOK_PATH.relative_to(ROOT)),
        "static": static_report,
        "runtime": runtime_report,
        "outputs": output_report,
    }
    AUDIT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
