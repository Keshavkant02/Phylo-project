from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import tempfile
import textwrap
import urllib.request
import zipfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

import h5py
import nbformat as nbf
import numpy as np


ROOT = Path(__file__).resolve().parent
SOURCE_DIR = ROOT / "tmp" / "atacama_qiime2_source"
CACHE_DIR = ROOT / "soil_16s_class_cache"
NOTEBOOK_PATH = ROOT / "soil_microbiome_16s_class_safe_colab.ipynb"

METADATA_PATH = SOURCE_DIR / "sample_metadata.tsv"
TABLE_QZA = SOURCE_DIR / "atacama-table.qza"
REP_SEQS_QZA = SOURCE_DIR / "atacama-rep-seqs.qza"

SOURCE_URLS = {
    "table": [
        # docs.qiime2.org mirrors are reachable in more locked-down environments
        # than data.qiime2.org, while serving the same tutorial artifact bytes here.
        "https://docs.qiime2.org/2024.10/data/tutorials/chimera/atacama-table.qza",
        "https://data.qiime2.org/2024.10/tutorials/chimera/atacama-table.qza",
    ],
    "rep_seqs": [
        "https://docs.qiime2.org/2024.10/data/tutorials/chimera/atacama-rep-seqs.qza",
        "https://data.qiime2.org/2024.10/tutorials/chimera/atacama-rep-seqs.qza",
    ],
    "metadata": [
        # The Atacama tutorial metadata currently lives on data.qiime2.org.
        # The builder keeps a local copy once fetched.
        "https://data.qiime2.org/2024.10/tutorials/atacama-soils/sample_metadata.tsv",
    ],
}

SOURCE_CANDIDATES = {
    "table": [TABLE_QZA, SOURCE_DIR / "docs_atacama_table.qza"],
    "rep_seqs": [REP_SEQS_QZA, SOURCE_DIR / "docs_atacama_rep_seqs.qza"],
    "metadata": [METADATA_PATH],
}

PREFIX = "goal2_atacama"


def dedent(text: str) -> str:
    return textwrap.dedent(text).strip() + "\n"


def read_qza_file(qza_path: Path, suffix: str) -> bytes:
    with zipfile.ZipFile(qza_path) as zf:
        matches = [name for name in zf.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f"Expected one {suffix} in {qza_path}, found {matches}")
        return zf.read(matches[0])


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def qza_has_member(path: Path, suffix: str) -> bool:
    if not path.exists() or path.stat().st_size < 100 or not zipfile.is_zipfile(path):
        return False
    try:
        with zipfile.ZipFile(path) as zf:
            return any(name.endswith(suffix) for name in zf.namelist())
    except zipfile.BadZipFile:
        return False


def metadata_has_required_columns(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 100:
        return False
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = set(reader.fieldnames or [])
        required = {"sample-id", "average-soil-relative-humidity", "vegetation", "percentcover"}
        return required.issubset(fields)
    except Exception:
        return False


def download_source(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".download")
    with urllib.request.urlopen(url, timeout=120) as response:
        tmp.write_bytes(response.read())
    tmp.replace(destination)


def resolve_source_file(kind: str, validator, destination: Path) -> tuple[Path, str]:
    for candidate in SOURCE_CANDIDATES[kind]:
        if validator(candidate):
            return candidate, "local"

    errors = []
    for url in SOURCE_URLS[kind]:
        try:
            download_source(url, destination)
            if validator(destination):
                return destination, url
            errors.append(f"{url} downloaded but failed validation")
        except Exception as exc:
            errors.append(f"{url}: {exc}")

    searched = ", ".join(str(path) for path in SOURCE_CANDIDATES[kind])
    raise FileNotFoundError(
        f"Could not resolve a real QIIME 2 source artifact for {kind}. "
        f"Searched local files: {searched}. Download attempts: {' | '.join(errors)}. "
        "Place the real artifact in tmp/atacama_qiime2_source and rerun; this builder does not synthesize counts."
    )


def resolve_source_artifacts() -> tuple[Path, Path, Path, dict[str, dict[str, str]]]:
    table_qza, table_source = resolve_source_file("table", lambda p: qza_has_member(p, "feature-table.biom"), TABLE_QZA)
    rep_qza, rep_source = resolve_source_file("rep_seqs", lambda p: qza_has_member(p, "dna-sequences.fasta"), REP_SEQS_QZA)
    metadata_path, metadata_source = resolve_source_file("metadata", metadata_has_required_columns, METADATA_PATH)
    resolution = {
        "table": {"path": str(table_qza.relative_to(ROOT)), "source": table_source, "sha256": sha256_file(table_qza)},
        "rep_seqs": {"path": str(rep_qza.relative_to(ROOT)), "source": rep_source, "sha256": sha256_file(rep_qza)},
        "metadata": {"path": str(metadata_path.relative_to(ROOT)), "source": metadata_source, "sha256": sha256_file(metadata_path)},
    }
    return metadata_path, table_qza, rep_qza, resolution


def read_biom_table(qza_path: Path) -> tuple[list[str], list[str], np.ndarray]:
    payload = read_qza_file(qza_path, "feature-table.biom")
    with tempfile.NamedTemporaryFile(delete=False, suffix=".biom") as tmp:
        tmp.write(payload)
        tmp_path = Path(tmp.name)
    try:
        with h5py.File(tmp_path, "r") as biom:
            feature_ids = [value.decode("utf-8") for value in biom["observation/ids"][:]]
            sample_ids = [value.decode("utf-8") for value in biom["sample/ids"][:]]
            data = biom["observation/matrix/data"][:]
            indices = biom["observation/matrix/indices"][:]
            indptr = biom["observation/matrix/indptr"][:]

        table = np.zeros((len(feature_ids), len(sample_ids)), dtype=float)
        for feature_index in range(len(feature_ids)):
            start = int(indptr[feature_index])
            end = int(indptr[feature_index + 1])
            table[feature_index, indices[start:end]] = data[start:end]
        return feature_ids, sample_ids, table
    finally:
        tmp_path.unlink(missing_ok=True)


def read_fasta_from_qza(qza_path: Path) -> dict[str, str]:
    text = read_qza_file(qza_path, "dna-sequences.fasta").decode("utf-8")
    records: dict[str, list[str]] = OrderedDict()
    current = None
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = []
        elif current is not None:
            records[current].append(line.strip().upper())
    return {key: "".join(chunks) for key, chunks in records.items()}


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: dict[str, dict[str, str]] = OrderedDict()
        for row in reader:
            sample_id = row["sample-id"]
            if sample_id == "#q2:types":
                continue
            rows[sample_id] = row
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def safe_float(value: object) -> float:
    try:
        if value in ("", None):
            return math.nan
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def shannon(counts: np.ndarray) -> float:
    total = float(np.sum(counts))
    if total <= 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(-(p * np.log(p)).sum())


def clean_taxonomy_rank(value: str, *, rank: str) -> str:
    value = value.strip()
    if "__" in value:
        value = value.split("__", 1)[1].strip()
    if not value or value.lower() in {"uncultured", "unidentified", "unknown", "unassigned"}:
        return ""
    if rank == "genus" and (value.endswith("aceae") or value.endswith("ales")):
        # QIIME/SILVA can repeat a family/order placeholder at genus level.
        return ""
    return value


def parse_taxonomy_string(taxonomy: str) -> dict[str, str]:
    ranks = {"phylum": "", "family": "", "genus": ""}
    for chunk in taxonomy.split(";"):
        value = chunk.strip()
        if value.startswith("p__"):
            ranks["phylum"] = clean_taxonomy_rank(value, rank="phylum")
        elif value.startswith("f__"):
            ranks["family"] = clean_taxonomy_rank(value, rank="family")
        elif value.startswith("g__"):
            ranks["genus"] = clean_taxonomy_rank(value, rank="genus")
    return ranks


def read_taxonomy_if_available() -> tuple[dict[str, dict[str, str]], str]:
    candidates = [
        CACHE_DIR / f"{PREFIX}_qiime_taxonomy.tsv",
        SOURCE_DIR / "taxonomy.tsv",
        SOURCE_DIR / "atacama-taxonomy.tsv",
        SOURCE_DIR / "atacama-taxonomy.qza",
        SOURCE_DIR / "taxonomy.qza",
        SOURCE_DIR / "docs_atacama_taxonomy.qza",
    ]
    for path in candidates:
        if not path.exists() or path.stat().st_size < 100:
            continue
        try:
            if path.suffix == ".tsv":
                text = path.read_text(encoding="utf-8")
            else:
                if not zipfile.is_zipfile(path):
                    continue
                with zipfile.ZipFile(path) as zf:
                    matches = [name for name in zf.namelist() if name.endswith("taxonomy.tsv")]
                    if len(matches) != 1:
                        continue
                    text = zf.read(matches[0]).decode("utf-8")
            reader = csv.DictReader(io.StringIO(text), delimiter="\t")
            taxonomy: dict[str, dict[str, str]] = {}
            for row in reader:
                feature_id = row.get("Feature ID") or row.get("feature-id") or row.get("FeatureID") or row.get("id")
                tax = row.get("Taxon") or row.get("taxonomy") or ""
                if not feature_id:
                    continue
                ranks = parse_taxonomy_string(tax)
                match = ranks["genus"] or ranks["family"] or "Unassigned at genus level"
                taxonomy[feature_id] = {
                    "phylum": ranks["phylum"] or "Unassigned",
                    "family": ranks["family"] or "",
                    "genus": ranks["genus"] or "",
                    "closest_taxonomic_match": match,
                }
            if taxonomy:
                return taxonomy, f"Taxonomy loaded from {path.name}."
        except Exception:
            continue

    static_assignment_path = CACHE_DIR / f"{PREFIX}_silva_static_taxonomy_assignments.csv"
    if static_assignment_path.exists() and static_assignment_path.stat().st_size > 100:
        try:
            with static_assignment_path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                taxonomy: dict[str, dict[str, str]] = {}
                for row in reader:
                    feature_id = row.get("qiime_feature_id") or ""
                    if not feature_id:
                        continue
                    taxonomy[feature_id] = {
                        "phylum": row.get("phylum") or "Unassigned",
                        "family": row.get("family") or "",
                        "genus": row.get("genus") or "",
                        "closest_taxonomic_match": row.get("closest_taxonomic_match") or "Unassigned at genus level",
                    }
            if taxonomy:
                return (
                    taxonomy,
                    "Fallback closest SILVA 138 515F/806R reference matches loaded from the local taxonomy cache.",
                )
        except Exception:
            pass
    return {}, "No valid SILVA taxonomy artifact was found locally; taxonomy is not inferred."


def format_float(value: float, digits: int = 5) -> float:
    if not math.isfinite(value):
        return math.nan
    return round(float(value), digits)


@dataclass(frozen=True)
class CachePaths:
    metadata: Path
    counts_top50: Path
    relative_top20: Path
    alpha: Path
    feature_key: Path
    sequences: Path
    manifest: Path
    readme: Path


def build_cache() -> CachePaths:
    metadata_path, table_qza, rep_seqs_qza, source_resolution = resolve_source_artifacts()

    CACHE_DIR.mkdir(exist_ok=True)
    metadata = read_metadata(metadata_path)
    feature_ids, sample_ids, table = read_biom_table(table_qza)
    sequences_by_feature = read_fasta_from_qza(rep_seqs_qza)
    taxonomy_by_feature, taxonomy_note = read_taxonomy_if_available()

    if len(sample_ids) != 61:
        raise ValueError(f"Goal 2 expects 61 samples; found {len(sample_ids)}")
    missing_metadata = sorted(set(sample_ids) - set(metadata))
    if missing_metadata:
        raise ValueError(f"Feature table samples missing metadata: {missing_metadata[:5]}")

    sample_totals = table.sum(axis=0)
    relative = np.divide(table, sample_totals, out=np.zeros_like(table), where=sample_totals > 0) * 100.0
    total_reads = table.sum(axis=1)
    prevalence = (table > 0).sum(axis=1)
    mean_relative = relative.mean(axis=1)
    max_relative = relative.max(axis=1)

    present_threshold = math.ceil(0.10 * len(sample_ids))
    eligible_for_q = [idx for idx, value in enumerate(prevalence) if value >= present_threshold]
    # The 61-sample source table has only nine ASVs present in >=10% of samples.
    # To keep the BH multiple-testing lesson honest and at the requested scale, use
    # the top 50 ASVs by prevalence and document that low-prevalence ASVs have lower power.
    q_indices = sorted(range(len(feature_ids)), key=lambda idx: (-prevalence[idx], -mean_relative[idx], feature_ids[idx]))[:50]
    top20_indices = list(np.argsort(-mean_relative)[:20])
    top12_indices = list(np.argsort(-mean_relative)[:12])
    top8_indices = list(np.argsort(-mean_relative)[:8])
    union_indices = sorted(set(q_indices) | set(top20_indices) | set(top12_indices) | set(top8_indices), key=lambda idx: (-mean_relative[idx], feature_ids[idx]))

    label_by_feature = {feature_ids[idx]: f"Atacama_ASV_{rank:02d}" for rank, idx in enumerate(union_indices, start=1)}
    feature_by_label = {label: feature_id for feature_id, label in label_by_feature.items()}

    metadata_rows: list[dict[str, object]] = []
    for sample_id in sample_ids:
        row = metadata[sample_id]
        metadata_rows.append(
            {
                "sample_id": sample_id,
                "transect_name": row["transect-name"],
                "site_name": row["site-name"],
                "depth": row["depth"],
                "elevation": row["elevation"],
                "average_soil_relative_humidity": row["average-soil-relative-humidity"],
                "vegetation": row["vegetation"],
                "percentcover": row["percentcover"],
            }
        )

    metadata_path = CACHE_DIR / f"{PREFIX}_sample_metadata.csv"
    write_csv(
        metadata_path,
        metadata_rows,
        [
            "sample_id",
            "transect_name",
            "site_name",
            "depth",
            "elevation",
            "average_soil_relative_humidity",
            "vegetation",
            "percentcover",
        ],
    )

    q_columns = [label_by_feature[feature_ids[idx]] for idx in q_indices]
    q_rows: list[dict[str, object]] = []
    for sample_pos, sample_id in enumerate(sample_ids):
        row: dict[str, object] = {"sample_id": sample_id}
        for idx in q_indices:
            row[label_by_feature[feature_ids[idx]]] = int(table[idx, sample_pos])
        row["total_reads"] = int(sample_totals[sample_pos])
        q_rows.append(row)
    counts_top50_path = CACHE_DIR / f"{PREFIX}_counts_top50.csv"
    write_csv(counts_top50_path, q_rows, ["sample_id", *q_columns, "total_reads"])

    top20_labels = [label_by_feature[feature_ids[idx]] for idx in top20_indices]
    relative_rows: list[dict[str, object]] = []
    for sample_pos, sample_id in enumerate(sample_ids):
        row = {"sample_id": sample_id}
        top_sum = 0.0
        for idx in top20_indices:
            label = label_by_feature[feature_ids[idx]]
            value = float(relative[idx, sample_pos])
            row[label] = format_float(value)
            top_sum += value
        row["Other"] = format_float(max(0.0, 100.0 - top_sum))
        row["total_reads"] = int(sample_totals[sample_pos])
        relative_rows.append(row)
    relative_top20_path = CACHE_DIR / f"{PREFIX}_relative_abundance_top20.csv"
    write_csv(relative_top20_path, relative_rows, ["sample_id", *top20_labels, "Other", "total_reads"])

    alpha_rows: list[dict[str, object]] = []
    for sample_pos, sample_id in enumerate(sample_ids):
        counts = table[:, sample_pos]
        alpha_rows.append(
            {
                "sample_id": sample_id,
                "total_reads": int(sample_totals[sample_pos]),
                "observed_asvs": int(np.sum(counts > 0)),
                "shannon_diversity": format_float(shannon(counts)),
            }
        )
    alpha_path = CACHE_DIR / f"{PREFIX}_alpha_diversity.csv"
    write_csv(alpha_path, alpha_rows, ["sample_id", "total_reads", "observed_asvs", "shannon_diversity"])

    q_feature_set = {feature_ids[idx] for idx in q_indices}
    top20_set = {feature_ids[idx] for idx in top20_indices}
    top12_set = {feature_ids[idx] for idx in top12_indices}
    top8_set = {feature_ids[idx] for idx in top8_indices}

    feature_rows: list[dict[str, object]] = []
    for idx in union_indices:
        feature_id = feature_ids[idx]
        label = label_by_feature[feature_id]
        tax = taxonomy_by_feature.get(
            feature_id,
            {
                "phylum": "Unassigned",
                "family": "",
                "genus": "",
                "closest_taxonomic_match": "Unassigned at genus level",
            },
        )
        feature_rows.append(
            {
                "asv": label,
                "qiime_feature_id": feature_id,
                "total_reads": int(total_reads[idx]),
                "prevalence_samples": int(prevalence[idx]),
                "mean_relative_abundance_percent": format_float(mean_relative[idx]),
                "max_relative_abundance_percent": format_float(max_relative[idx]),
                "sequence_length": len(sequences_by_feature.get(feature_id, "")),
                "closest_taxonomic_match": tax["closest_taxonomic_match"] or "Unassigned at genus level",
                "phylum": tax["phylum"] or "Unassigned",
                "family": tax["family"],
                "genus": tax["genus"],
                "in_q_value_top50": feature_id in q_feature_set,
                "in_abundance_top20": feature_id in top20_set,
                "in_tree_top12": feature_id in top12_set,
                "in_alignment_top8": feature_id in top8_set,
            }
        )
    feature_key_path = CACHE_DIR / f"{PREFIX}_feature_key.csv"
    write_csv(
        feature_key_path,
        feature_rows,
        [
            "asv",
            "qiime_feature_id",
            "total_reads",
            "prevalence_samples",
            "mean_relative_abundance_percent",
            "max_relative_abundance_percent",
            "sequence_length",
            "closest_taxonomic_match",
            "phylum",
            "family",
            "genus",
            "in_q_value_top50",
            "in_abundance_top20",
            "in_tree_top12",
            "in_alignment_top8",
        ],
    )

    sequences_path = CACHE_DIR / f"{PREFIX}_rep_seqs_top50_union.fasta"
    with sequences_path.open("w", encoding="utf-8", newline="\n") as handle:
        for label, feature_id in sorted(feature_by_label.items()):
            sequence = sequences_by_feature.get(feature_id, "")
            if not sequence:
                raise ValueError(f"Missing representative sequence for {feature_id}")
            handle.write(f">{label} qiime_feature_id={feature_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    manifest = {
        "title": "Goal 2 Atacama soil ASV cache",
        "source_table": "https://data.qiime2.org/2024.10/tutorials/chimera/atacama-table.qza",
        "source_rep_seqs": "https://data.qiime2.org/2024.10/tutorials/chimera/atacama-rep-seqs.qza",
        "source_metadata": "https://data.qiime2.org/2024.10/tutorials/atacama-soils/sample_metadata.tsv",
        "source_context": "QIIME 2 2024.10 Atacama soil tutorial and q2-vsearch chimera tutorial artifacts.",
        "source_resolution": source_resolution,
        "data_mode": "real_qiime2_artifacts_only_no_synthetic_counts",
        "taxonomy_note": taxonomy_note,
        "sample_count": len(sample_ids),
        "feature_count_full_table": len(feature_ids),
        "present_threshold_samples": present_threshold,
        "asvs_at_or_above_10_percent_prevalence": len(eligible_for_q),
        "q_value_asvs": len(q_indices),
        "abundance_asvs": len(top20_indices),
        "tree_asvs": len(top12_indices),
        "alignment_asvs": len(top8_indices),
        "scientific_note": "Taxonomy labels come from a real QIIME/SILVA taxonomy artifact when available, with the nearest-reference cache used only as a fallback. Labels are closest matches, not species proof.",
    }
    manifest_path = CACHE_DIR / f"{PREFIX}_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    readme_path = CACHE_DIR / "ATACAMA_GOAL2_README.md"
    readme_path.write_text(
        dedent(
            """
            # Goal 2 Atacama Soil ASV Cache

            This cache supports the Atacama-only student Colab in `soil_microbiome_16s_class_safe_colab.ipynb`.
            It is derived from QIIME 2 2024.10 Atacama tutorial artifacts and contains no synthetic ASVs.

            Source files used locally:

            - a validated QIIME 2 feature table artifact containing `feature-table.biom`
            - a validated QIIME 2 representative-sequence artifact containing `dna-sequences.fasta`
            - Atacama sample metadata containing `average-soil-relative-humidity`

            The builder can use already-local files in `tmp/atacama_qiime2_source/`, or fetch the
            official QIIME tutorial artifacts when the network is available. Broken 404 files and
            non-zip `.qza` placeholders are rejected. If real artifacts cannot be found, the builder
            stops instead of synthesizing counts.

            The notebook uses fixed subsets for readability:

            - top 50 ASVs by prevalence for association tests
            - top 20 ASVs by mean relative abundance for abundance plots
            - top 12 ASVs by mean relative abundance for the UPGMA tree
            - top 8 ASVs by mean relative abundance for the alignment heatmap

            Prevalence policy: the source artifact has only nine ASVs present in at least 10% of samples,
            so the notebook uses the top 50 by prevalence for the BH correction lesson and treats the
            lowest-prevalence ASVs cautiously.

            Taxonomy policy: the preferred student-facing cache reads
            `goal2_atacama_qiime_taxonomy.tsv`, produced by QIIME 2 `feature-classifier classify-sklearn`
            with the SILVA 138 Naive Bayes classifier. The builder can also read local QIIME taxonomy
            artifacts from `tmp/atacama_qiime2_source/`. The older
            `goal2_atacama_silva_static_taxonomy_assignments.csv` nearest-reference cache is kept only as a
            documented fallback. Taxonomy remains a closest-match label, not species proof.
            """
        ),
        encoding="utf-8",
    )

    return CachePaths(
        metadata=metadata_path,
        counts_top50=counts_top50_path,
        relative_top20=relative_top20_path,
        alpha=alpha_path,
        feature_key=feature_key_path,
        sequences=sequences_path,
        manifest=manifest_path,
        readme=readme_path,
    )


def parse_fasta_text(text: str) -> OrderedDict[str, str]:
    records: OrderedDict[str, list[str]] = OrderedDict()
    current = None
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = []
        elif current is not None:
            records[current].append(line.strip())
    return OrderedDict((key, "".join(value)) for key, value in records.items())


def make_setup_cell(cache_files: dict[str, str]) -> str:
    return dedent(
        f'''
        import io
        import json
        import math
        import random
        import warnings
        from collections import Counter
        from dataclasses import dataclass
        from pathlib import Path

        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        from IPython.display import HTML, Markdown, display
        from matplotlib import patches
        from scipy import stats
        try:
            from statsmodels.stats.multitest import multipletests
        except ModuleNotFoundError:
            # Fallback only if a runtime lacks statsmodels. The intended path is statsmodels' fdr_bh.
            def multipletests(pvals, method="fdr_bh"):
                pvals = np.asarray(pvals, dtype=float)
                order = np.argsort(pvals)
                adjusted = np.empty_like(pvals)
                running = 1.0
                n = len(pvals)
                for rank_from_end, idx in enumerate(order[::-1], start=1):
                    rank = n - rank_from_end + 1
                    running = min(running, pvals[idx] * n / rank)
                    adjusted[idx] = min(running, 1.0)
                return adjusted < 0.05, adjusted, None, None

        warnings.filterwarnings("ignore", category=RuntimeWarning)

        CACHE_FILES = {cache_files!r}

        OKABE_ITO = {{
            "blue": "#0072B2",
            "orange": "#E69F00",
            "green": "#009E73",
            "pink": "#CC79A7",
            "sky": "#56B4E9",
            "vermillion": "#D55E00",
            "yellow": "#F0E442",
            "black": "#000000",
        }}
        BASE_COLORS = {{"A": "#56B4E9", "C": "#009E73", "G": "#E69F00", "T": "#CC79A7", "-": "#EEEEEE", "N": "#EEEEEE"}}
        SAMPLE_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000"]

        plt.rcParams.update({{
            "figure.dpi": 130,
            "savefig.dpi": 180,
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.frameon": False,
        }})

        def cache_text(name):
            local = Path("soil_16s_class_cache") / name
            if local.exists():
                return local.read_text(encoding="utf-8")
            return CACHE_FILES[name]

        def read_cache_csv(name):
            return pd.read_csv(io.StringIO(cache_text(name)))

        def parse_fasta(text):
            records = {{}}
            current = None
            chunks = []
            for line in text.splitlines():
                if not line.strip():
                    continue
                if line.startswith(">"):
                    if current is not None:
                        records[current] = "".join(chunks)
                    current = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line.strip().upper())
            if current is not None:
                records[current] = "".join(chunks)
            return records

        def add_caption(fig, text):
            fig.text(0.01, -0.04, text, ha="left", va="top", fontsize=9, style="italic", color="#4D4D4D")

        def clean_axes(ax):
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            return ax

        def styled_table(df, width_px=880):
            return (
                df.style.hide(axis="index")
                .set_table_styles([
                    {{"selector": "table", "props": [("border-collapse", "collapse"), ("width", f"{{width_px}}px"), ("font-size", "13px")]}},
                    {{"selector": "th", "props": [("text-align", "left"), ("border-bottom", "1px solid #999"), ("padding", "6px 8px")]}},
                    {{"selector": "td", "props": [("padding", "6px 8px"), ("border-bottom", "1px solid #e6e6e6")]}},
                    {{"selector": "tbody tr:nth-child(odd)", "props": [("background-color", "#fafafa")]}},
                ])
            )

        def fmt_p(value):
            if pd.isna(value):
                return ""
            value = float(value)
            return f"{{value:.2e}}" if value < 0.001 else f"{{value:.3f}}"

        def short_label(asv):
            return asv.replace("_", " ")

        def sequence_distance(seq_a, seq_b):
            compared = 0
            differences = 0
            for a, b in zip(seq_a, seq_b):
                if a not in "ACGT" or b not in "ACGT":
                    continue
                compared += 1
                differences += int(a != b)
            return differences / compared if compared else math.nan

        metadata = read_cache_csv("goal2_atacama_sample_metadata.csv")
        counts_top50 = read_cache_csv("goal2_atacama_counts_top50.csv")
        relative_top20 = read_cache_csv("goal2_atacama_relative_abundance_top20.csv")
        alpha_diversity = read_cache_csv("goal2_atacama_alpha_diversity.csv")
        feature_key = read_cache_csv("goal2_atacama_feature_key.csv")
        manifest = json.loads(cache_text("goal2_atacama_manifest.json"))
        sequences = parse_fasta(cache_text("goal2_atacama_rep_seqs_top50_union.fasta"))

        for col in ["average_soil_relative_humidity", "percentcover"]:
            metadata[col] = pd.to_numeric(metadata[col], errors="coerce")
        for col in ["mean_relative_abundance_percent", "max_relative_abundance_percent", "prevalence_samples"]:
            feature_key[col] = pd.to_numeric(feature_key[col], errors="coerce")

        feature_key = feature_key.sort_values("mean_relative_abundance_percent", ascending=False).reset_index(drop=True)
        top20_asvs = feature_key.query("in_abundance_top20 == True")["asv"].tolist()
        top12_asvs = feature_key.query("in_tree_top12 == True")["asv"].tolist()
        top8_asvs = feature_key.query("in_alignment_top8 == True")["asv"].tolist()
        q_asvs = feature_key.query("in_q_value_top50 == True")["asv"].tolist()

        # ASV cascade - chosen for readability at each step:
        #   top 50 by prevalence  -> q-value tests (statistical power)
        #   top 20 by mean abundance -> abundance plots (readable bars)
        #   top 12 by mean abundance -> UPGMA tree (readable tip labels)
        #   top 8  by mean abundance -> alignment heatmap (readable bases)
        # The 61-sample source table has fewer than 50 ASVs at >=10% prevalence,
        # so the q-value section keeps the requested top-50 scale and interprets low-prevalence results cautiously.

        asv_color = {{asv: SAMPLE_COLORS[i % len(SAMPLE_COLORS)] for i, asv in enumerate(top20_asvs)}}
        asv_color["Other"] = "#BDBDBD"
        '''
    )


def make_tree_helpers_cell() -> str:
    return dedent(
        '''
        @dataclass
        class TreeNode:
            name: str | None = None
            left: object | None = None
            right: object | None = None
            left_length: float = 0.0
            right_length: float = 0.0

            @property
            def is_leaf(self):
                return self.name is not None

        def mammal_tree(rotated=False):
            dog_wolf = TreeNode(left=TreeNode("dog"), right=TreeNode("wolf"), left_length=1.0, right_length=1.0)
            dog_wolf_fox = TreeNode(left=dog_wolf, right=TreeNode("fox"), left_length=1.0, right_length=2.0)
            cat_lion = TreeNode(left=TreeNode("cat"), right=TreeNode("lion"), left_length=2.0, right_length=2.0)
            bear_cat_lion = TreeNode(left=TreeNode("bear"), right=cat_lion, left_length=3.0, right_length=1.5)
            if rotated:
                dog_wolf = TreeNode(left=TreeNode("wolf"), right=TreeNode("dog"), left_length=1.0, right_length=1.0)
                dog_wolf_fox = TreeNode(left=TreeNode("fox"), right=dog_wolf, left_length=2.0, right_length=1.0)
                cat_lion = TreeNode(left=TreeNode("lion"), right=TreeNode("cat"), left_length=2.0, right_length=2.0)
                bear_cat_lion = TreeNode(left=cat_lion, right=TreeNode("bear"), left_length=1.5, right_length=3.0)
            return TreeNode(left=dog_wolf_fox, right=bear_cat_lion, left_length=1.0, right_length=0.5)

        def assign_y(node, leaf_order):
            lookup = {name: i for i, name in enumerate(leaf_order)}
            ys = {}
            def walk(current):
                if current.is_leaf:
                    ys[id(current)] = lookup[current.name]
                    return ys[id(current)]
                left_y = walk(current.left)
                right_y = walk(current.right)
                ys[id(current)] = (left_y + right_y) / 2
                return ys[id(current)]
            walk(node)
            return ys

        def draw_tree(ax, node, leaf_order, label_map=None, x0=0.0, y_lookup=None, color="#333333", label_offset=0.08):
            if y_lookup is None:
                y_lookup = assign_y(node, leaf_order)
            def walk(current, x):
                y = y_lookup[id(current)]
                if current.is_leaf:
                    text = label_map.get(current.name, current.name) if label_map else current.name
                    ax.text(x + label_offset, y, text, va="center", ha="left", fontsize=9)
                    return
                children = [(current.left, current.left_length), (current.right, current.right_length)]
                child_ys = []
                for child, length in children:
                    cy = y_lookup[id(child)]
                    cx = x + length
                    ax.plot([x, cx], [cy, cy], color=color, lw=1.4)
                    child_ys.append(cy)
                    walk(child, cx)
                ax.plot([x, x], [min(child_ys), max(child_ys)], color=color, lw=1.4)
            walk(node, x0)
            ax.set_ylim(-0.6, len(leaf_order) - 0.4)
            ax.set_yticks([])
            ax.set_xlabel("branch length (DNA-difference units)")
            clean_axes(ax)
            ax.invert_yaxis()
            ax.margins(x=0.12)
            return ax
        '''
    )


def code_cell_tree_layouts() -> str:
    return dedent(
        '''
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharex=True)
        orders = [
            ["dog", "wolf", "fox", "bear", "cat", "lion"],
            ["fox", "dog", "wolf", "bear", "cat", "lion"],
            ["cat", "lion", "bear", "fox", "wolf", "dog"],
        ]
        titles = ["Rectangular layout", "Ladderized layout", "Alternative tip order"]
        for ax, order, title in zip(axes, orders, titles):
            draw_tree(ax, mammal_tree(rotated=False), order)
            ax.set_title(title, loc="left")
            ax.text(0.02, 1.02, "Which two tips are closest relatives?", transform=ax.transAxes, ha="left", va="bottom", fontsize=9, color="#4D4D4D")
        fig.suptitle("The same mammal tree can look different without changing relationships.", x=0.01, ha="left", fontsize=12)
        add_caption(fig, "Answer: dog and wolf are sister taxa in all three layouts; reading left-to-right across tip order is not enough.")
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_rotation_diagram() -> str:
    return dedent(
        '''
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharex=True)
        draw_tree(axes[0], mammal_tree(rotated=False), ["dog", "wolf", "fox", "bear", "cat", "lion"])
        axes[0].set_title("Before node rotation", loc="left")
        draw_tree(axes[1], mammal_tree(rotated=True), ["fox", "wolf", "dog", "lion", "cat", "bear"])
        axes[1].set_title("After rotating internal nodes", loc="left")
        for ax in axes:
            ax.scatter([1.0, 3.5], [0.5, 4.5], s=42, facecolor="white", edgecolor=OKABE_ITO["green"], zorder=3)
            ax.text(0.03, 0.02, "Tips = living organisms\\nInternal nodes = inferred common ancestors", transform=ax.transAxes, ha="left", va="bottom", fontsize=9, color="#4D4D4D")
        axes[0].annotate("MRCA of dog and wolf", xy=(1.0, 0.5), xytext=(1.7, 1.6), arrowprops=dict(arrowstyle="-", color="#666666"), fontsize=9)
        axes[0].annotate("Sister groups", xy=(3.5, 4.5), xytext=(3.7, 3.6), arrowprops=dict(arrowstyle="-", color="#666666"), fontsize=9)
        fig.suptitle("Rotating a node changes the drawing, not the ancestry hypothesis.", x=0.01, ha="left", fontsize=12)
        add_caption(fig, "Branch length is read from the horizontal axis here; in this notebook it means DNA difference, not time.")
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_dataset_summary() -> str:
    return dedent(
        '''
        summary = pd.DataFrame([
            {"What": "Samples", "Value": f"{len(metadata)}", "Why it matters": "Each row is one soil sample."},
            {"What": "ASVs available for tests", "Value": f"{len(q_asvs)}", "Why it matters": "Enough ASVs to ask which patterns track humidity or vegetation."},
            {"What": "Tree ASVs", "Value": f"{len(top12_asvs)}", "Why it matters": "Small enough for readable tip labels."},
            {"What": "Metadata used", "Value": "humidity, vegetation", "Why it matters": "These describe the soil environment."},
            {"What": "Taxonomic names", "Value": manifest["taxonomy_note"], "Why it matters": "Names are closest-reference labels, not proof of exact species."},
        ])
        display(styled_table(summary, width_px=940))
        '''
    )


def sequence_distance(seq_a: str, seq_b: str) -> float:
    compared = 0
    differences = 0
    for a, b in zip(seq_a, seq_b):
        if a not in "ACGT" or b not in "ACGT":
            continue
        compared += 1
        differences += int(a != b)
    return differences / compared if compared else math.nan


def code_cell_alignment() -> str:
    return dedent(
        '''
        alignment_asvs = top8_asvs
        matrix = np.array([[base for base in sequences[asv]] for asv in alignment_asvs])

        def identity_fraction(column):
            bases = [base for base in column if base in "ACGT"]
            if not bases:
                return 1.0
            counts = Counter(bases)
            return max(counts.values()) / len(bases)

        window_width = 60
        best_start = 0
        best_score = -1
        for start in range(0, matrix.shape[1] - window_width + 1):
            window = matrix[:, start:start + window_width]
            identity = np.array([identity_fraction(window[:, col]) for col in range(window_width)])
            variable = np.sum(identity < 0.875)
            conserved = np.sum(identity >= 0.875)
            gap_penalty = np.sum(window == "-")
            score = variable + 0.12 * conserved - 0.2 * gap_penalty
            if variable > 0 and score > best_score:
                best_start = start
                best_score = score

        window = matrix[:, best_start:best_start + window_width]
        identity = np.array([identity_fraction(window[:, col]) for col in range(window_width)])
        alpha_by_column = np.where(identity >= 0.875, 0.30, 1.0)

        fig_height = 0.5 * len(alignment_asvs) + 1.2
        fig, ax = plt.subplots(figsize=(12, fig_height))
        for row_idx, asv in enumerate(alignment_asvs):
            for col_idx, base in enumerate(window[row_idx]):
                color = BASE_COLORS.get(base, "#EEEEEE")
                ax.add_patch(patches.Rectangle((col_idx, row_idx), 1, 1, facecolor=color, edgecolor="none", alpha=float(alpha_by_column[col_idx])))

        labels = [short_label(asv) for asv in alignment_asvs]
        ax.set_xlim(0, window_width)
        ax.set_ylim(0, len(alignment_asvs))
        ax.set_yticks(np.arange(len(alignment_asvs)) + 0.5)
        ax.set_yticklabels(labels, fontfamily="DejaVu Sans Mono", fontsize=10)
        for label in ax.get_yticklabels():
            label.set_horizontalalignment("right")
        ax.tick_params(axis="y", pad=18, length=0)
        tick_positions = [0.5, window_width / 2, window_width - 0.5]
        tick_labels = [best_start + 1, best_start + window_width // 2, best_start + window_width]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        ax.tick_params(axis="x", top=False, bottom=True, length=3)
        ax.set_xlabel("aligned marker-window column")
        ax.set_title("Aligned 16S marker window - variable columns carry the phylogenetic signal.", loc="left")
        ax.invert_yaxis()
        for spine in ax.spines.values():
            spine.set_visible(False)
        add_caption(fig, "Colored cells = DNA bases (A blue, C green, G orange, T pink). Faded columns = conserved across all sequences; bright columns = where the sequences differ.")
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_distance_matrix() -> str:
    return dedent(
        '''
        distance_asvs = top12_asvs
        distance_matrix = pd.DataFrame(index=distance_asvs, columns=distance_asvs, dtype=float)
        for a in distance_asvs:
            for b in distance_asvs:
                distance_matrix.loc[a, b] = sequence_distance(sequences[a], sequences[b])

        fig, ax = plt.subplots(figsize=(7.2, 6.1))
        shown = distance_matrix.loc[distance_asvs, distance_asvs]
        im = ax.imshow(shown.values, cmap="viridis", vmin=0, vmax=np.nanmax(shown.values))
        ax.set_xticks(range(len(distance_asvs)))
        ax.set_yticks(range(len(distance_asvs)))
        ax.set_xticklabels([short_label(x).replace("Atacama ", "") for x in distance_asvs], rotation=45, ha="right")
        ax.set_yticklabels([short_label(x).replace("Atacama ", "") for x in distance_asvs])
        if len(distance_asvs) <= 12:
            for i in range(len(distance_asvs)):
                for j in range(len(distance_asvs)):
                    ax.text(j, i, f"{shown.values[i, j]:.2f}", ha="center", va="center", fontsize=7, color="white" if shown.values[i, j] > 0.18 else "#222222")
        ax.set_title("Pairwise DNA distance among top Atacama ASVs.", loc="left")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("sequence-distance units")
        clean_axes(ax)
        add_caption(fig, "Small distances mark ASV sequences with more similar DNA letters in comparable columns.")
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_upgma_tree() -> str:
    return dedent(
        '''
        @dataclass
        class UPGMANode:
            name: str | None = None
            left: object | None = None
            right: object | None = None
            height: float = 0.0
            left_length: float = 0.0
            right_length: float = 0.0
            members: tuple[str, ...] = ()

            @property
            def is_leaf(self):
                return self.name is not None

        def build_upgma(distance_df):
            clusters = {name: UPGMANode(name=name, height=0.0, members=(name,)) for name in distance_df.index}
            next_id = 1
            while len(clusters) > 1:
                keys = list(clusters)
                best = None
                best_distance = math.inf
                for i, a in enumerate(keys):
                    for b in keys[i + 1:]:
                        vals = [distance_df.loc[x, y] for x in clusters[a].members for y in clusters[b].members]
                        d = float(np.mean(vals))
                        if d < best_distance:
                            best = (a, b)
                            best_distance = d
                a, b = best
                left = clusters.pop(a)
                right = clusters.pop(b)
                height = best_distance / 2.0
                node = UPGMANode(
                    left=left,
                    right=right,
                    height=height,
                    left_length=max(0.0, height - left.height),
                    right_length=max(0.0, height - right.height),
                    members=tuple(sorted(left.members + right.members)),
                )
                clusters[f"node_{next_id}"] = node
                next_id += 1
            return next(iter(clusters.values()))

        def leaves_in_order(node):
            if node.is_leaf:
                return [node.name]
            return leaves_in_order(node.left) + leaves_in_order(node.right)

        def newick(node):
            if node.is_leaf:
                return node.name
            return f"({newick(node.left)}:{node.left_length:.6f},{newick(node.right)}:{node.right_length:.6f})"

        upgma_tree = build_upgma(distance_matrix)
        tree_newick = newick(upgma_tree) + ";"
        leaf_order = leaves_in_order(upgma_tree)
        y_lookup = {name: i for i, name in enumerate(leaf_order)}

        def node_y(node):
            if node.is_leaf:
                return y_lookup[node.name]
            return (node_y(node.left) + node_y(node.right)) / 2

        tax_lookup = feature_key.set_index("asv")["closest_taxonomic_match"].to_dict()
        fig_height = max(4.8, 0.48 * len(leaf_order) + 1.5)
        fig, ax = plt.subplots(figsize=(9.5, fig_height))

        def draw_upgma(node, x):
            y = node_y(node)
            if node.is_leaf:
                match = tax_lookup.get(node.name, "Unassigned at genus level")
                suffix = "unassigned" if match == "Unassigned at genus level" else match
                ax.text(x + 0.01, y, f"{short_label(node.name)} ({suffix})", va="center", ha="left", fontsize=9)
                return
            children = [(node.left, node.left_length), (node.right, node.right_length)]
            child_ys = []
            for child, length in children:
                cy = node_y(child)
                cx = x + length
                ax.plot([x, cx], [cy, cy], color="#333333", lw=1.3)
                child_ys.append(cy)
                draw_upgma(child, cx)
            ax.plot([x, x], [min(child_ys), max(child_ys)], color="#333333", lw=1.3)

        draw_upgma(upgma_tree, 0.0)
        ax.set_ylim(-0.6, len(leaf_order) - 0.4)
        ax.invert_yaxis()
        ax.set_yticks([])
        ax.set_xlabel("branch length (sequence-distance units)")
        ax.set_title("UPGMA tree for the top Atacama ASVs.", loc="left")
        clean_axes(ax)
        ax.margins(x=0.25)
        ax.set_xlim(left=0)
        add_caption(fig, "Tips are ASVs. Short paths between tips suggest closer sequence relatedness, not exact species identity.")
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_abundance_table() -> str:
    return dedent(
        '''
        abundance_table = (
            feature_key.query("in_abundance_top20 == True")
            .loc[:, ["asv", "closest_taxonomic_match", "mean_relative_abundance_percent", "max_relative_abundance_percent", "prevalence_samples"]]
            .rename(columns={
                "asv": "ASV",
                "closest_taxonomic_match": "Closest taxonomic match",
                "mean_relative_abundance_percent": "Mean relative abundance (%)",
                "max_relative_abundance_percent": "Maximum relative abundance (%)",
                "prevalence_samples": "Samples detected (of 61)",
            })
        )
        abundance_table["Mean relative abundance (%)"] = abundance_table["Mean relative abundance (%)"].round(2)
        abundance_table["Maximum relative abundance (%)"] = abundance_table["Maximum relative abundance (%)"].round(2)
        display(styled_table(abundance_table, width_px=1000))
        '''
    )


def code_cell_abundance_plot() -> str:
    return dedent(
        '''
        plot_df = relative_top20.merge(metadata[["sample_id", "average_soil_relative_humidity", "vegetation"]], on="sample_id")
        plot_df = plot_df.sort_values("average_soil_relative_humidity").reset_index(drop=True)
        stack_cols = [col for col in top20_asvs if col in plot_df.columns] + ["Other"]

        fig, ax = plt.subplots(figsize=(12, 4.6))
        bottom = np.zeros(len(plot_df))
        x = np.arange(len(plot_df))
        for col in stack_cols:
            values = plot_df[col].to_numpy(dtype=float)
            ax.bar(x, values, bottom=bottom, width=0.88, color=asv_color.get(col, "#BDBDBD"), edgecolor="white", linewidth=0.15, label=short_label(col))
            bottom += values

        ax.set_ylim(0, 100)
        ax.set_xlim(-0.5, len(plot_df) - 0.5)
        ax.set_ylabel("relative abundance (%)")
        ax.set_xlabel("samples ordered from drier to wetter soil")
        ax.set_xticks([])
        ax.set_title("Top ASVs across Atacama soil samples ordered by humidity.", loc="left")
        ax.text(0, -9, "drier", ha="left", va="top", fontsize=9, color="#4D4D4D")
        ax.text(len(plot_df) - 1, -9, "wetter", ha="right", va="top", fontsize=9, color="#4D4D4D")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:10], labels[:10], bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, title="first 10 ASVs")
        clean_axes(ax)
        add_caption(fig, "Each vertical bar is one soil sample; color shows which ASVs make up the sample after the top 20 are separated from Other.")
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_alpha_table() -> str:
    return dedent(
        '''
        alpha_plot = alpha_diversity.merge(metadata[["sample_id", "average_soil_relative_humidity", "vegetation"]], on="sample_id")
        alpha_summary = (
            alpha_plot.groupby("vegetation", as_index=False)
            .agg(
                Samples=("sample_id", "count"),
                **{
                    "Mean observed ASVs": ("observed_asvs", "mean"),
                    "Mean Shannon diversity": ("shannon_diversity", "mean"),
                    "Mean humidity (%)": ("average_soil_relative_humidity", "mean"),
                },
            )
        )
        alpha_summary["Mean observed ASVs"] = alpha_summary["Mean observed ASVs"].round(1)
        alpha_summary["Mean Shannon diversity"] = alpha_summary["Mean Shannon diversity"].round(2)
        alpha_summary["Mean humidity (%)"] = alpha_summary["Mean humidity (%)"].round(1)
        alpha_summary = alpha_summary.rename(columns={"vegetation": "Vegetation"})
        display(styled_table(alpha_summary, width_px=760))
        '''
    )


def code_cell_alpha_plot() -> str:
    return dedent(
        '''
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8), sharex=True)
        metrics = [("observed_asvs", "observed ASVs"), ("shannon_diversity", "Shannon diversity")]
        for ax, (metric, label) in zip(axes, metrics):
            x = alpha_plot["average_soil_relative_humidity"].to_numpy(dtype=float)
            y = alpha_plot[metric].to_numpy(dtype=float)
            ax.scatter(x, y, s=28, color=OKABE_ITO["blue"], alpha=0.78, edgecolor="white", linewidth=0.3)
            valid = np.isfinite(x) & np.isfinite(y)
            if valid.sum() >= 3:
                fit = np.polyfit(x[valid], y[valid], 1)
                xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
                ax.plot(xs, fit[0] * xs + fit[1], color="#666666", lw=1.2, alpha=0.7)
            ax.set_xlabel("average soil relative humidity (%)")
            ax.set_ylabel(label)
            ax.set_title(label, loc="left")
            clean_axes(ax)
        fig.suptitle("Alpha diversity across the humidity gradient.", x=0.01, ha="left", fontsize=12)
        add_caption(fig, "Observed ASVs count types; Shannon diversity increases when a sample has many types with more even abundances.")
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_association_plot() -> str:
    return dedent(
        '''
        count_cols = [col for col in counts_top50.columns if col.startswith("Atacama_ASV_")]
        q_counts = counts_top50[["sample_id", *count_cols]].merge(metadata[["sample_id", "average_soil_relative_humidity", "vegetation"]], on="sample_id")
        raw_counts = q_counts[count_cols].astype(float)
        # CLR uses a 0.5 pseudo-count so zero counts can be logged without creating infinite values.
        logged = np.log(raw_counts + 0.5)
        clr = logged.sub(logged.mean(axis=1), axis=0)
        sample_totals = counts_top50.set_index("sample_id")["total_reads"].reindex(q_counts["sample_id"]).to_numpy(dtype=float)
        rel = raw_counts.div(sample_totals, axis=0) * 100.0

        humidity = q_counts["average_soil_relative_humidity"].astype(float)
        vegetation = q_counts["vegetation"].astype(str).str.lower()
        rows = []
        for asv in count_cols:
            rho, p_humidity = stats.spearmanr(clr[asv], humidity, nan_policy="omit")
            yes_values = clr.loc[vegetation == "yes", asv]
            no_values = clr.loc[vegetation == "no", asv]
            if len(yes_values) > 0 and len(no_values) > 0:
                u_stat, p_veg = stats.mannwhitneyu(yes_values, no_values, alternative="two-sided")
            else:
                p_veg = np.nan
            mean_yes = rel.loc[vegetation == "yes", asv].mean()
            mean_no = rel.loc[vegetation == "no", asv].mean()
            log2_fc = float(np.log2((mean_yes + 0.001) / (mean_no + 0.001)))
            mean_abundance = rel[asv].mean()
            rows.append({"ASV": asv, "Variable": "Humidity", "Effect": float(rho), "Raw p-value": float(p_humidity), "Mean abundance (%)": float(mean_abundance)})
            rows.append({"ASV": asv, "Variable": "Vegetation", "Effect": log2_fc, "Raw p-value": float(p_veg), "Mean abundance (%)": float(mean_abundance)})

        association_results = pd.DataFrame(rows)
        q_values = []
        for variable, group in association_results.groupby("Variable", sort=False):
            _, adjusted, _, _ = multipletests(group["Raw p-value"].fillna(1.0).to_numpy(dtype=float), method="fdr_bh")
            q_values.extend(pd.Series(adjusted, index=group.index).items())
        q_lookup = dict(q_values)
        association_results["BH q-value"] = association_results.index.map(q_lookup).astype(float)
        association_results = association_results.merge(
            feature_key[["asv", "closest_taxonomic_match", "phylum", "genus", "prevalence_samples"]].rename(columns={"asv": "ASV"}),
            on="ASV",
            how="left",
        )
        association_results["Significant?"] = np.where(association_results["BH q-value"] < 0.05, "yes", "-")

        humidity_results = association_results.query("Variable == 'Humidity'").copy()
        fig, axes = plt.subplots(1, 2, figsize=(4, 3), sharey=True)
        rng = np.random.default_rng(7)
        for ax, col, title in zip(axes, ["Raw p-value", "BH q-value"], ["Raw p-values", "BH q-values"]):
            y = humidity_results[col].to_numpy(dtype=float)
            x = rng.normal(0, 0.035, size=len(y))
            ax.scatter(x, y, s=18, color=OKABE_ITO["blue"], alpha=0.62, edgecolor="none")
            ax.axhline(0.05, color="#666666", ls="--", lw=1)
            ax.set_xlim(-0.18, 0.18)
            ax.set_xticks([])
            ax.set_title(title, loc="left", fontsize=10)
            clean_axes(ax)
        axes[0].set_ylabel("value")
        fig.suptitle("Raw p-values compared with BH q-values.", x=0.02, ha="left", fontsize=11)
        add_caption(fig, "Multiple-testing correction is conservative on purpose - it is the cost of asking many questions at once.")
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_q_table() -> str:
    return dedent(
        '''
        table_rows = []
        for variable in ["Humidity", "Vegetation"]:
            subset = association_results.query("Variable == @variable").sort_values("BH q-value").head(10).copy()
            for _, row in subset.iterrows():
                if variable == "Humidity":
                    effect_label = f"{row['Effect']:.2f}"
                    effect_col = "Spearman rho"
                else:
                    arrow = "higher with vegetation" if row["Effect"] > 0 else "higher without vegetation"
                    effect_label = f"{row['Effect']:.2f} ({arrow})"
                    effect_col = "Effect size"
                table_rows.append({
                    "Metadata variable": variable,
                    "ASV": row["ASV"],
                    "Closest taxonomic match": row["closest_taxonomic_match"],
                    "Mean abundance (%)": round(row["Mean abundance (%)"], 2),
                    "Samples detected (of 61)": int(row["prevalence_samples"]),
                    "Spearman rho": effect_label if variable == "Humidity" else "",
                    "Effect size": effect_label if variable == "Vegetation" else "",
                    "Raw p-value": fmt_p(row["Raw p-value"]),
                    "BH q-value": fmt_p(row["BH q-value"]),
                    "Significant?": row["Significant?"],
                })

        q_table = pd.DataFrame(table_rows)

        def significant_style(row):
            if row["Significant?"] == "yes":
                return ["border-left: 3px solid #009E73"] + [""] * (len(row) - 1)
            return [""] * len(row)

        def sig_color(value):
            if value == "yes":
                return "color: #009E73; font-weight: 700"
            if value == "-":
                return "color: #999999"
            return ""

        def q_background(value):
            try:
                numeric = float(str(value).replace("e", "E"))
            except ValueError:
                return ""
            strength = 1.0 - min(max(numeric / 0.05, 0.0), 1.0)
            return f"background-color: rgba(0, 158, 115, {0.06 + 0.18 * strength:.3f})"

        display(
            q_table.style.hide(axis="index")
            .apply(significant_style, axis=1)
            .map(sig_color, subset=["Significant?"])
            .map(q_background, subset=["BH q-value"])
            .set_table_styles([
                {"selector": "table", "props": [("border-collapse", "collapse"), ("width", "1080px"), ("font-size", "12px")]},
                {"selector": "th", "props": [("text-align", "left"), ("border-bottom", "1px solid #999"), ("padding", "6px 7px")]},
                {"selector": "td", "props": [("padding", "6px 7px"), ("border-bottom", "1px solid #e6e6e6")]},
                {"selector": "tbody tr:nth-child(odd)", "props": [("background-color", "#fafafa")]},
            ])
        )
        '''
    )


def code_cell_lollipop() -> str:
    return dedent(
        '''
        # The full top-50 test table includes low-prevalence ASVs, but the teaching plot
        # only labels discoveries that meet the original >=10% prevalence threshold.
        min_prevalence_for_plot = int(manifest["present_threshold_samples"])
        significant = association_results.query("`BH q-value` < 0.05 and prevalence_samples >= @min_prevalence_for_plot").copy()
        panels = ["Humidity", "Vegetation"]
        fig, axes = plt.subplots(1, 2, figsize=(11, max(3.4, 0.38 * max(1, significant.groupby("Variable").size().max() if not significant.empty else 1) + 1.4)))
        for ax, variable in zip(axes, panels):
            subset = significant.query("Variable == @variable").copy()
            if subset.empty:
                ax.text(0.5, 0.5, "No ASVs below q < 0.05", ha="center", va="center", transform=ax.transAxes, color="#666666")
                ax.set_axis_off()
                continue
            subset = subset.sort_values("Effect")
            y = np.arange(len(subset))
            colors = [OKABE_ITO["blue"] if row["phylum"] in ("", "Unassigned") else OKABE_ITO["green"] for _, row in subset.iterrows()]
            sizes = 45 + 18 * np.sqrt(subset["Mean abundance (%)"].clip(lower=0.01))
            ax.axvline(0, color="#777777", lw=1)
            ax.hlines(y, 0, subset["Effect"], color="#777777", lw=1.1)
            ax.scatter(subset["Effect"], y, s=sizes, color=colors, alpha=0.88, edgecolor="white", linewidth=0.4, zorder=3)
            labels = []
            for _, row in subset.iterrows():
                match = row["genus"] if isinstance(row["genus"], str) and row["genus"] else row.get("closest_taxonomic_match", "unassigned")
                if not isinstance(match, str) or not match:
                    match = "unassigned"
                labels.append(f"{short_label(row['ASV'])} ({match})")
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=8)
            x_right = max(subset["Effect"].max(), 0) + 0.08
            for yi, (_, row) in enumerate(subset.iterrows()):
                ax.text(x_right, yi, f"q={fmt_p(row['BH q-value'])}", ha="left", va="center", fontsize=8, color="#999999")
            ax.set_xlabel("Spearman rho" if variable == "Humidity" else "log2 fold-change")
            ax.set_title(variable, loc="left")
            clean_axes(ax)
            ax.margins(x=0.22)
        fig.suptitle("ASVs with BH q-values below 0.05 and enough prevalence to interpret.", x=0.01, y=0.98, ha="left", fontsize=12)
        add_caption(fig, "Each bar is one ASV. Humidity bars to the right increase with wetter soil; vegetation bars to the right are higher with vegetation. Dot size = overall abundance. Only ASVs with BH q < 0.05 and >=10% prevalence shown.")
        plt.tight_layout(rect=[0, 0.06, 1, 0.91])
        plt.show()
        '''
    )


def build_notebook(cache_paths: CachePaths) -> None:
    cache_files = {
        cache_paths.metadata.name: cache_paths.metadata.read_text(encoding="utf-8"),
        cache_paths.counts_top50.name: cache_paths.counts_top50.read_text(encoding="utf-8"),
        cache_paths.relative_top20.name: cache_paths.relative_top20.read_text(encoding="utf-8"),
        cache_paths.alpha.name: cache_paths.alpha.read_text(encoding="utf-8"),
        cache_paths.feature_key.name: cache_paths.feature_key.read_text(encoding="utf-8"),
        cache_paths.sequences.name: cache_paths.sequences.read_text(encoding="utf-8"),
        cache_paths.manifest.name: cache_paths.manifest.read_text(encoding="utf-8"),
    }

    nb = nbf.v4.new_notebook()
    nb.metadata.update(
        {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
            "colab": {"name": NOTEBOOK_PATH.name, "provenance": []},
        }
    )

    cells = []
    md = nbf.v4.new_markdown_cell
    code = nbf.v4.new_code_cell

    cells.append(
        md(
            dedent(
                """
                # Soil 16S phylogeny: Atacama ASVs

                ## Section 1: Story hook - Atacama Desert

                The Atacama Desert is often described as the driest non-polar desert on Earth. Some places there have recorded zero rainfall across decades, yet microbial life persists in the soil.

                Today we ask one question from three angles: **Who lives there, where are they abundant, and how are they related?**
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 2: Tree-thinking intro (mammals only)

                A phylogenetic tree is a hypothesis about relatedness. Before looking at microbes, use familiar mammals to practice reading trees: dog, wolf, fox, bear, cat, and lion.

                First question: if the same tree is drawn three ways, do the closest relatives change?
                """
            )
        )
    )
    cells.append(code(make_setup_cell(cache_files) + "\n" + make_tree_helpers_cell()))
    cells.append(code(code_cell_tree_layouts()))
    cells.append(
        md(
            dedent(
                """
                In all three drawings, dog and wolf remain sister taxa because they share the most recent common ancestor with each other. Tip order alone can trick your eye; the branching pattern is what carries the evidence.

                Think: If a tree is rotated around an internal node, what should stay the same?
                """
            )
        )
    )
    cells.append(code(code_cell_rotation_diagram()))
    cells.append(
        md(
            dedent(
                """
                *A tip is a present-day organism or sequence. An internal node is an inferred common ancestor, not a sample whose DNA we directly measured. A sister group is a pair of lineages that share an immediate common ancestor; the MRCA is the most recent common ancestor for the taxa you are tracing.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 3: Atacama dataset story

                The rest of the notebook uses real 16S amplicon data from the QIIME 2 Atacama soil tutorial and its 2024.10 Atacama teaching artifacts. These are soil samples collected across a humidity and aridity gradient, with metadata recording whether vegetation was present near each sample.
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 4: What is an ASV?

                An **ASV** is a precise DNA sequence pattern found after cleaning 16S sequencing reads. An ASV is not automatically a species.

                Keep three ideas separate: **abundance** means how much of an ASV is in a sample; **taxonomic match** means what known group the sequence resembles; **evolutionary relatedness** means how sequences cluster in a tree.
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 5: Load cached Atacama data

                The notebook starts from prepared Atacama tables so every student sees the same data when pressing Run all. The student-facing summary below shows the scale of the dataset without dumping raw metadata.
                """
            )
        )
    )
    cells.append(code(code_cell_dataset_summary()))
    cells.append(
        md(
            dedent(
                """
                *The dataset has enough samples to compare dry and wetter soils, but the tree will use a smaller ASV subset so the labels stay readable.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 6: Alignment of representative ASV sequences

                Alignment means putting DNA letters into comparable columns. Conserved columns are mostly the same across ASVs; variable columns are where sequence differences become visible.
                """
            )
        )
    )
    cells.append(code(code_cell_alignment()))
    cells.append(
        md(
            dedent(
                """
                *ASVs that look similar across the bright columns are more closely related. The phylogenetic tree (section 8) formalizes this intuition.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 7: Distance matrix

                Small distance means more similar DNA sequence. Here, each number compares two ASVs in the same aligned 16S marker region.
                """
            )
        )
    )
    cells.append(code(code_cell_distance_matrix()))
    cells.append(
        md(
            dedent(
                """
                *The darkest cells mark the smallest distances; those pairs are the first candidates for close sequence relatedness.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 8: UPGMA tree (the payoff)

                Tips are ASVs. Branch length shows sequence difference. Branches that meet recently, with a short path between them, suggest closer sequence relatedness, which is our best evidence of evolutionary relatedness here but not proof of exact species identity.
                """
            )
        )
    )
    cells.append(code(code_cell_upgma_tree()))
    cells.append(
        md(
            dedent(
                """
                *This UPGMA tree turns the distance matrix into a readable hypothesis: nearby tips have more similar 16S sequences than tips far apart on the tree.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 9: Relative abundance - what is actually in these samples?

                The tree asks how sequences are related. Abundance asks a different question: which ASVs make up each soil sample, and does that composition change from drier to wetter soils?
                """
            )
        )
    )
    cells.append(code(code_cell_abundance_table()))
    cells.append(code(code_cell_abundance_plot()))
    cells.append(
        md(
            dedent(
                """
                *A very abundant ASV is not automatically the most evolutionarily unusual one; abundance and relatedness answer different questions.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 10: Alpha diversity

                Alpha diversity asks how diverse one sample is: how many ASV types are present, and how evenly distributed they are. We use observed ASVs for richness and Shannon diversity for richness plus evenness.
                """
            )
        )
    )
    cells.append(code(code_cell_alpha_table()))
    cells.append(code(code_cell_alpha_plot()))
    cells.append(
        md(
            dedent(
                """
                *Use the trend lines as a visual guide, then describe the direction carefully: humid samples may be richer or more even, but the plot does not identify exact species.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 11: BH-corrected association tests

                We tested 50 ASVs to see if abundance changes with humidity. With a p < 0.05 cutoff, we would expect about 50 x 0.05 = 2.5 false positives by random chance alone, even if nothing is truly associated. Benjamini-Hochberg (BH) correction adjusts for this. A q-value of 0.05 means we expect about 5% of discoveries below that threshold to be false alarms.
                """
            )
        )
    )
    cells.append(code(code_cell_association_plot()))
    cells.append(code(code_cell_q_table()))
    cells.append(code(code_cell_lollipop()))
    cells.append(
        md(
            dedent(
                """
                *BH q-values summarize abundance-versus-metadata tests. They are not tree branch support values, and they do not say an ASV is a proved species.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 12: Final student report

                Fill in the report using the section numbers named in each question.

                1. Which 3 ASVs are most abundant in your samples? Refer to section 9 table.  
                   [your answer]

                2. Which ASVs are significantly associated with humidity? Refer to section 11 lollipop.  
                   [your answer]

                3. Which ASVs are significantly associated with vegetation? Refer to section 11 lollipop.  
                   [your answer]

                4. What does alpha diversity suggest about humid versus arid samples? Refer to section 10.  
                   [your answer]

                5. In the UPGMA tree, which ASVs cluster closest together? Refer to section 8.  
                   [your answer]

                6. Are closely related ASVs from section 8 also similar in abundance from section 9 or in humidity association from section 11?  
                   [your answer]

                7. Synthesis: An ASV can be (a) very abundant, (b) statistically associated with humidity, and (c) closely related to another ASV in the tree - and these are three different things. Explain in 2-3 sentences why these are three different ideas and why all three matter.  
                   [your answer]
                """
            )
        )
    )

    nb["cells"] = cells
    if not (28 <= len(cells) <= 36):
        raise AssertionError(f"Goal 2 requires 28-36 cells; built {len(cells)}")
    NOTEBOOK_PATH.write_text(nbf.writes(nb), encoding="utf-8")


def main() -> None:
    cache_paths = build_cache()
    build_notebook(cache_paths)
    print(f"Wrote {NOTEBOOK_PATH.name}")
    print(f"Wrote Goal 2 cache files under {CACHE_DIR}")


if __name__ == "__main__":
    main()
