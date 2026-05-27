from __future__ import annotations

import base64
import csv
import hashlib
import io
import json
import math
import re
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
DARWIN_INTRO_IMAGE_PATH = CACHE_DIR / "darwin_tree_1837_intro.png"


def dedent(text: str) -> str:
    return textwrap.dedent(text).strip() + "\n"


def image_data_uri(path: Path, mime_type: str = "image/png") -> str:
    if not path.exists() or path.stat().st_size < 100:
        raise FileNotFoundError(f"Missing cached image for notebook embedding: {path}")
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime_type};base64,{encoded}"


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
    lowered = value.lower()
    blocked = {"uncultured", "unidentified", "unknown", "unassigned", "incertae_sedis", "incertae sedis"}
    if not value or lowered in blocked:
        return ""
    if "uncultured" in lowered or "unassigned" in lowered:
        return ""
    if rank == "genus" and (value.endswith("aceae") or value.endswith("ales")):
        # QIIME/SILVA can repeat a family/order placeholder at genus level.
        return ""
    if re.match(r"^[A-Za-z]*\d+[A-Za-z0-9-]*$", value) or re.match(r"^\d", value):
        # Codes such as wb1-P19 or 67-14 are honest taxonomy strings, but not useful names.
        return ""
    return value


def parse_taxonomy_string(taxonomy: str) -> dict[str, str]:
    ranks = {"phylum": "", "class": "", "order": "", "family": "", "genus": ""}
    for chunk in taxonomy.split(";"):
        value = chunk.strip()
        if value.startswith("p__"):
            ranks["phylum"] = clean_taxonomy_rank(value, rank="phylum")
        elif value.startswith("c__"):
            ranks["class"] = clean_taxonomy_rank(value, rank="class")
        elif value.startswith("o__"):
            ranks["order"] = clean_taxonomy_rank(value, rank="order")
        elif value.startswith("f__"):
            ranks["family"] = clean_taxonomy_rank(value, rank="family")
        elif value.startswith("g__"):
            ranks["genus"] = clean_taxonomy_rank(value, rank="genus")
    return ranks


def taxonomy_display_label(ranks: dict[str, str]) -> str:
    if ranks.get("genus"):
        return ranks["genus"]
    for rank in ("family", "order", "class", "phylum"):
        if ranks.get(rank):
            return f"{ranks[rank]} ({rank})"
    return "Unassigned"


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
                match = taxonomy_display_label(ranks)
                taxonomy[feature_id] = {
                    "phylum": ranks["phylum"] or "Unassigned",
                    "class": ranks["class"],
                    "order": ranks["order"],
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
                        "class": row.get("class") or "",
                        "order": row.get("order") or "",
                        "family": row.get("family") or "",
                        "genus": row.get("genus") or "",
                        "closest_taxonomic_match": row.get("closest_taxonomic_match") or "Unassigned",
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
    counts_retained: Path
    relative_top20: Path
    alpha: Path
    feature_key: Path
    sequences_retained: Path
    manifest: Path
    readme: Path


def build_cache() -> CachePaths:
    metadata_path, table_qza, rep_seqs_qza, source_resolution = resolve_source_artifacts()

    CACHE_DIR.mkdir(exist_ok=True)
    metadata = read_metadata(metadata_path)
    feature_ids, raw_sample_ids, raw_table = read_biom_table(table_qza)
    sequences_by_feature = read_fasta_from_qza(rep_seqs_qza)
    taxonomy_by_feature, taxonomy_note = read_taxonomy_if_available()

    if len(raw_sample_ids) != 61:
        raise ValueError(f"Goal 2 expects 61 raw samples; found {len(raw_sample_ids)}")
    missing_metadata = sorted(set(raw_sample_ids) - set(metadata))
    if missing_metadata:
        raise ValueError(f"Feature table samples missing metadata: {missing_metadata[:5]}")

    raw_sample_totals = raw_table.sum(axis=0)
    qc_sample_indices: list[int] = []
    dropped_zero_reads = 0
    dropped_low_reads_nonzero = 0
    dropped_missing_metadata = 0
    for sample_pos, sample_id in enumerate(raw_sample_ids):
        row = metadata[sample_id]
        humidity = safe_float(row.get("average-soil-relative-humidity"))
        vegetation = (row.get("vegetation") or "").strip().lower()
        has_metadata = math.isfinite(humidity) and vegetation in {"yes", "no"}
        total = float(raw_sample_totals[sample_pos])
        if total < 100:
            if total == 0:
                dropped_zero_reads += 1
            else:
                dropped_low_reads_nonzero += 1
            continue
        if not has_metadata:
            dropped_missing_metadata += 1
            continue
        qc_sample_indices.append(sample_pos)

    sample_ids = [raw_sample_ids[idx] for idx in qc_sample_indices]
    table = raw_table[:, qc_sample_indices]
    sample_totals = table.sum(axis=0)
    relative = np.divide(table, sample_totals, out=np.zeros_like(table), where=sample_totals > 0) * 100.0
    total_reads = table.sum(axis=1)
    prevalence = (table > 0).sum(axis=1)
    mean_relative = relative.mean(axis=1)
    max_relative = relative.max(axis=1)

    if len(sample_ids) != 46:
        raise ValueError(f"Final brief expects 46 QC-passed samples; found {len(sample_ids)}")
    if (dropped_zero_reads, dropped_low_reads_nonzero, dropped_missing_metadata) != (5, 7, 3):
        raise ValueError(
            "Unexpected sample-QC drop counts: "
            f"zero={dropped_zero_reads}, low_nonzero={dropped_low_reads_nonzero}, missing_metadata={dropped_missing_metadata}"
        )

    prevalence_threshold = 3
    retained_indices = [idx for idx, value in enumerate(prevalence) if value >= prevalence_threshold]
    if len(retained_indices) != 37:
        raise ValueError(f"Final brief expects 37 ASVs present in >=3 QC samples; found {len(retained_indices)}")

    # ASV cascade - chosen for readability at each step:
    #   prevalence >= 3 samples first -> q-value tests (all retained ASVs)
    #   top 20 by mean relative abundance -> abundance plots (readable bars)
    #   top 12 by mean relative abundance -> UPGMA tree (readable tip labels)
    #   top 8  by mean relative abundance -> alignment heatmap (readable bases)
    # Mean relative abundance is used instead of total reads to avoid deep-sample bias.
    sorted_retained = sorted(retained_indices, key=lambda idx: (-mean_relative[idx], feature_ids[idx]))
    q_indices = sorted_retained
    top20_indices = sorted_retained[:20]
    top12_indices = sorted_retained[:12]
    top8_indices = sorted_retained[:8]
    union_indices = sorted_retained

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
                "total_reads": int(sample_totals[len(metadata_rows)]),
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
            "total_reads",
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
    counts_retained_path = CACHE_DIR / f"{PREFIX}_counts_retained_asvs.csv"
    write_csv(counts_retained_path, q_rows, ["sample_id", *q_columns, "total_reads"])

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
                "class": "",
                "order": "",
                "family": "",
                "genus": "",
                "closest_taxonomic_match": "Unassigned",
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
                "closest_taxonomic_match": tax["closest_taxonomic_match"] or "Unassigned",
                "phylum": tax["phylum"] or "Unassigned",
                "class": tax.get("class", ""),
                "order": tax.get("order", ""),
                "family": tax["family"],
                "genus": tax["genus"],
                "in_q_value_tests": feature_id in q_feature_set,
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
            "class",
            "order",
            "family",
            "genus",
            "in_q_value_tests",
            "in_abundance_top20",
            "in_tree_top12",
            "in_alignment_top8",
        ],
    )

    sequences_path = CACHE_DIR / f"{PREFIX}_rep_seqs_retained_asvs.fasta"
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
        "qiime2_version": "Amplicon 2024.10.1",
        "classifier": "SILVA 138 99% OTUs full-length Naive Bayes classifier",
        "classifier_sha256": "c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616",
        "raw_sample_count": len(raw_sample_ids),
        "sample_count": len(sample_ids),
        "samples_dropped_by_qc": len(raw_sample_ids) - len(sample_ids),
        "samples_dropped_zero_reads": dropped_zero_reads,
        "samples_dropped_low_reads_nonzero": dropped_low_reads_nonzero,
        "samples_dropped_missing_metadata": dropped_missing_metadata,
        "feature_count_full_table": len(feature_ids),
        "prevalence_filter_samples": prevalence_threshold,
        "retained_asvs": len(retained_indices),
        "q_value_asvs": len(q_indices),
        "abundance_asvs": len(top20_indices),
        "tree_asvs": len(top12_indices),
        "alignment_asvs": len(top8_indices),
        "scientific_note": "Taxonomy labels come from a real QIIME/SILVA taxonomy artifact when available. Labels are closest taxonomic matches, not species proof.",
    }
    manifest_path = CACHE_DIR / f"{PREFIX}_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")

    readme_path = CACHE_DIR / "ATACAMA_GOAL2_README.md"
    readme_path.write_text(
        dedent(
            """
            # Goal 2 Atacama Soil ASV Cache

            This cache supports the Atacama-only student Colab in `soil_microbiome_16s_class_safe_colab.ipynb`.
            It is derived from real QIIME 2 2024.10 Atacama tutorial artifacts and contains no synthetic ASVs.

            Source files used locally:

            - a validated QIIME 2 feature table artifact containing `feature-table.biom`
            - a validated QIIME 2 representative-sequence artifact containing `dna-sequences.fasta`
            - Atacama sample metadata containing `average-soil-relative-humidity`

            The builder can use already-local files in `tmp/atacama_qiime2_source/`, or fetch the
            official QIIME tutorial artifacts when the network is available. Broken 404 files and
            non-zip `.qza` placeholders are rejected. If real artifacts cannot be found, the builder
            stops instead of synthesizing counts.

            The notebook applies sample quality control first:

            - raw output: 401 ASVs across 61 samples
            - keep samples with at least 100 reads and complete humidity/vegetation metadata: 46 samples
            - keep ASVs present in at least three QC-passed samples: 37 ASVs

            The notebook then uses fixed subsets for readability:

            - all 37 retained ASVs for association tests
            - top 20 retained ASVs by mean relative abundance for abundance plots
            - top 12 retained ASVs by mean relative abundance for the UPGMA tree
            - top 8 retained ASVs by mean relative abundance for the alignment heatmap

            Sparseness policy: the 10% tutorial subsample is shallow. The notebook does not hide that;
            it makes quality control and conservative statistical claims part of the lesson.

            Taxonomy policy: the preferred student-facing cache reads
            `goal2_atacama_qiime_taxonomy.tsv`, produced by QIIME 2 `feature-classifier classify-sklearn`
            with the SILVA 138 Naive Bayes classifier. The builder can also read local QIIME taxonomy
            artifacts from `tmp/atacama_qiime2_source/`. Taxonomy remains a closest-match label, not species proof.
            """
        ),
        encoding="utf-8",
        newline="\n",
    )

    return CachePaths(
        metadata=metadata_path,
        counts_retained=counts_retained_path,
        relative_top20=relative_top20_path,
        alpha=alpha_path,
        feature_key=feature_key_path,
        sequences_retained=sequences_path,
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
        from matplotlib.lines import Line2D
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
            "figure.dpi": 110,
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

        def asv_tax_label(asv):
            match = tax_label_lookup.get(asv, "Unassigned")
            return f"{{short_label(asv)}} ({{match}})"

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
        counts_retained = read_cache_csv("goal2_atacama_counts_retained_asvs.csv")
        relative_top20 = read_cache_csv("goal2_atacama_relative_abundance_top20.csv")
        alpha_diversity = read_cache_csv("goal2_atacama_alpha_diversity.csv")
        feature_key = read_cache_csv("goal2_atacama_feature_key.csv")
        manifest = json.loads(cache_text("goal2_atacama_manifest.json"))
        sequences = parse_fasta(cache_text("goal2_atacama_rep_seqs_retained_asvs.fasta"))

        for col in ["average_soil_relative_humidity", "percentcover"]:
            metadata[col] = pd.to_numeric(metadata[col], errors="coerce")
        for col in ["mean_relative_abundance_percent", "max_relative_abundance_percent", "prevalence_samples"]:
            feature_key[col] = pd.to_numeric(feature_key[col], errors="coerce")

        feature_key = feature_key.sort_values("mean_relative_abundance_percent", ascending=False).reset_index(drop=True)
        top20_asvs = feature_key.query("in_abundance_top20 == True")["asv"].tolist()
        top12_asvs = feature_key.query("in_tree_top12 == True")["asv"].tolist()
        top8_asvs = feature_key.query("in_alignment_top8 == True")["asv"].tolist()
        q_asvs = feature_key.query("in_q_value_tests == True")["asv"].tolist()
        tax_label_lookup = feature_key.set_index("asv")["closest_taxonomic_match"].to_dict()
        phylum_lookup = feature_key.set_index("asv")["phylum"].fillna("Unassigned").to_dict()
        phyla = [p for p in feature_key["phylum"].fillna("Unassigned").unique().tolist()]
        phylum_color = {{phylum: SAMPLE_COLORS[i % len(SAMPLE_COLORS)] for i, phylum in enumerate(phyla)}}
        phylum_color["Unassigned"] = "#7F7F7F"

        # ASV cascade - chosen for readability at each step:
        #   sample QC first: >=100 reads and complete humidity/vegetation metadata
        #   ASV prevalence filter: present in >=3 QC-passed samples
        #   all 37 retained ASVs -> q-value tests
        #   top 20 retained by mean relative abundance -> abundance plots
        #   top 12 retained by mean relative abundance -> UPGMA tree
        #   top 8 retained by mean relative abundance -> alignment heatmap
        # CLR transform later uses a 0.5 pseudo-count before log transform.
        # Mean relative abundance is used instead of total reads to avoid deep-sample bias.

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

        def draw_tree(ax, node, leaf_order, label_map=None, x0=0.0, y_lookup=None, color="#333333", label_offset=0.08, xlabel="branch length units"):
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
            ax.set_xlabel(xlabel)
            clean_axes(ax)
            ax.invert_yaxis()
            ax.margins(x=0.12)
            return ax
        '''
    )


def code_cell_annotated_mammal_tree() -> str:
    return dedent(
        '''
        fig, ax = plt.subplots(figsize=(10.5, 5.0))
        label_map = {name: name.title() for name in ["dog", "wolf", "fox", "bear", "cat", "lion"]}
        draw_tree(ax, mammal_tree(rotated=False), ["dog", "wolf", "fox", "bear", "cat", "lion"], label_map=label_map, xlabel="branch length units")

        # Known coordinates from the teaching tree's branch lengths and leaf order.
        dog_tip = (3.0, 0.0)
        wolf_tip = (3.0, 1.0)
        fox_tip = (3.0, 2.0)
        dog_wolf_mrca = (2.0, 0.5)
        root = (0.0, 2.5)

        ax.scatter([dog_tip[0], wolf_tip[0], dog_wolf_mrca[0], root[0]], [dog_tip[1], wolf_tip[1], dog_wolf_mrca[1], root[1]],
                   s=[36, 36, 58, 58], facecolor="white", edgecolor=OKABE_ITO["green"], linewidth=1.4, zorder=4)
        ax.plot([3.28, 3.28], [0, 1], color=OKABE_ITO["green"], lw=1.1)
        ax.plot([3.21, 3.28], [0, 0], color=OKABE_ITO["green"], lw=1.1)
        ax.plot([3.21, 3.28], [1, 1], color=OKABE_ITO["green"], lw=1.1)

        ax.annotate("Tip: an organism being compared", xy=dog_tip, xytext=(3.35, -0.35),
                    arrowprops=dict(arrowstyle="-", color="#777777", lw=0.9), fontsize=9, color="#333333")
        ax.annotate("Branch: one line of descent", xy=(1.25, 1.25), xytext=(0.55, 0.35),
                    arrowprops=dict(arrowstyle="-", color="#777777", lw=0.9), fontsize=9, color="#333333")
        ax.annotate("Internal node: inferred common ancestor", xy=dog_wolf_mrca, xytext=(0.55, 1.75),
                    arrowprops=dict(arrowstyle="-", color="#777777", lw=0.9), fontsize=9, color="#333333")
        ax.annotate("Sister taxa share an immediate ancestor", xy=(3.28, 0.5), xytext=(3.55, 0.72),
                    arrowprops=dict(arrowstyle="-", color="#777777", lw=0.9), fontsize=9, color="#333333")
        ax.annotate("MRCA of all six mammals", xy=root, xytext=(0.35, 3.15),
                    arrowprops=dict(arrowstyle="-", color="#777777", lw=0.9), fontsize=9, color="#333333")
        ax.set_title("A phylogenetic tree is read by tracing branches back to common ancestors.", loc="left")
        add_caption(fig, "Relatedness is about common ancestry: dog and wolf share a more recent common ancestor with each other than either shares with fox.")
        plt.tight_layout()
        plt.show()
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
            draw_tree(ax, mammal_tree(rotated=False), order, xlabel="branch length units")
            ax.set_title(title, loc="left")
        fig.suptitle("The same mammal tree can look different without changing relationships.", x=0.01, ha="left", fontsize=12)
        add_caption(fig, "Answer: dog and wolf are sister taxa in all three layouts; reading left-to-right across tip order is not evidence of relatedness.")
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_rotation_diagram() -> str:
    return dedent(
        '''
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharex=True)
        draw_tree(axes[0], mammal_tree(rotated=False), ["dog", "wolf", "fox", "bear", "cat", "lion"], xlabel="branch length units")
        axes[0].set_title("Before node rotation", loc="left")
        draw_tree(axes[1], mammal_tree(rotated=True), ["fox", "wolf", "dog", "lion", "cat", "bear"], xlabel="branch length units")
        axes[1].set_title("After rotating internal nodes", loc="left")
        axes[0].scatter([2.0, 2.0], [0.5, 4.5], s=42, facecolor="white", edgecolor=OKABE_ITO["green"], zorder=3)
        axes[1].scatter([2.0, 2.0], [1.5, 3.5], s=42, facecolor="white", edgecolor=OKABE_ITO["green"], zorder=3)
        axes[0].annotate("MRCA of dog and wolf", xy=(2.0, 0.5), xytext=(2.25, 1.35), arrowprops=dict(arrowstyle="-", color="#666666"), fontsize=9)
        axes[0].annotate("Cat and lion are sister taxa", xy=(2.0, 4.5), xytext=(2.45, 3.8), arrowprops=dict(arrowstyle="-", color="#666666"), fontsize=9)
        fig.suptitle("Rotating a node changes the drawing, not the ancestry hypothesis.", x=0.01, ha="left", fontsize=12)
        add_caption(fig, "Rotation changes the drawing, not the ancestry hypothesis. In the microbial tree later, branch length means DNA sequence difference.")
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_dataset_summary() -> str:
    return dedent(
        '''
        humidity_min = metadata["average_soil_relative_humidity"].min()
        humidity_max = metadata["average_soil_relative_humidity"].max()
        vegetation_counts = metadata["vegetation"].str.lower().value_counts().to_dict()
        transects = ", ".join(sorted(metadata["transect_name"].dropna().unique()))
        summary = pd.DataFrame([
            {"Measure": "Samples (post-QC)", "Value": f"{len(metadata)}", "Meaning": "Each row is one usable soil sample."},
            {"Measure": "Samples dropped by QC", "Value": f"{manifest['samples_dropped_by_qc']} (5 zero reads, 7 more with <100 reads, 3 missing metadata)", "Meaning": "Very shallow or incomplete samples are unreliable for ASV counts."},
            {"Measure": "ASVs (post-prevalence filter)", "Value": f"{len(q_asvs)}", "Meaning": "These ASVs appeared in at least 3 QC-passed samples."},
            {"Measure": "Sequencing region", "Value": "16S V4 (~252 bp)", "Meaning": "A short marker region, not a whole genome."},
            {"Measure": "Humidity range (post-QC)", "Value": f"{humidity_min:.2f}-{humidity_max:.2f}%", "Meaning": "The samples still span dry to very humid soil."},
            {"Measure": "Vegetation groups", "Value": f"yes={vegetation_counts.get('yes', 0)}, no={vegetation_counts.get('no', 0)}", "Meaning": "Two metadata groups used later for association tests."},
            {"Measure": "Transects", "Value": transects, "Meaning": "The samples come from both Atacama transects in this tutorial subset."},
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

        labels = [asv_tax_label(asv) for asv in alignment_asvs]
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
        ax.set_title("Aligned 16S marker window — variable columns carry the phylogenetic signal.", loc="left")
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

        fig, ax = plt.subplots(figsize=(8.6, 6.8))
        shown = distance_matrix.loc[distance_asvs, distance_asvs]
        im = ax.imshow(shown.values, cmap="viridis", vmin=0, vmax=np.nanmax(shown.values))
        ax.set_xticks(range(len(distance_asvs)))
        ax.set_yticks(range(len(distance_asvs)))
        labels = [asv_tax_label(x).replace("Atacama ", "") for x in distance_asvs]
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels(labels, fontsize=7)
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

        fig_height = max(5.2, 0.52 * len(leaf_order) + 1.7)
        fig, ax = plt.subplots(figsize=(14, fig_height))

        def draw_upgma(node, x):
            y = node_y(node)
            if node.is_leaf:
                color = phylum_color.get(phylum_lookup.get(node.name, "Unassigned"), "#7F7F7F")
                ax.scatter([x], [y], s=32, color=color, edgecolor="white", linewidth=0.4, zorder=3)
                ax.text(x + 0.012, y, asv_tax_label(node.name), va="center", ha="left", fontsize=9)
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
        used_phyla = []
        for asv in leaf_order:
            phylum = phylum_lookup.get(asv, "Unassigned")
            if phylum not in used_phyla:
                used_phyla.append(phylum)
        handles = [
            Line2D([0], [0], marker="o", color="none", markerfacecolor=phylum_color.get(phylum, "#7F7F7F"),
                   markeredgecolor="white", markersize=7, label=phylum)
            for phylum in used_phyla
        ]
        ax.legend(handles=handles, title="Phylum", bbox_to_anchor=(1.01, 0.5), loc="center left", fontsize=8, title_fontsize=9)
        ax.set_xlim(0, max(upgma_tree.height * 2.9, upgma_tree.height + 0.08))
        add_caption(fig, "Tips are ASVs. Short paths between tips suggest closer sequence relatedness in the V4 marker region, not exact species identity.")
        plt.tight_layout(rect=[0, 0, 0.82, 1])
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
                "prevalence_samples": "Samples detected (of 46)",
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
        legend_pairs = list(zip(handles[:10], labels[:10]))
        if "Other" in stack_cols:
            other_index = stack_cols.index("Other")
            legend_pairs.append((handles[other_index], "Other"))
        ax.legend([h for h, _ in legend_pairs], [l for _, l in legend_pairs], bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, title="top 10 + Other")
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
            colors = alpha_plot["vegetation"].str.lower().map({"no": OKABE_ITO["orange"], "yes": OKABE_ITO["green"]}).fillna(OKABE_ITO["blue"])
            ax.scatter(x, y, s=30, color=colors, alpha=0.82, edgecolor="white", linewidth=0.35)
            valid = np.isfinite(x) & np.isfinite(y)
            if valid.sum() >= 3:
                fit = np.polyfit(x[valid], y[valid], 1)
                xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
                ax.plot(xs, fit[0] * xs + fit[1], color="#666666", lw=1.1, ls="--", alpha=0.65)
                rho, p_value = stats.spearmanr(x[valid], y[valid])
                ax.text(0.03, 0.95, f"Spearman rho={rho:.2f}\\np={fmt_p(p_value)}", transform=ax.transAxes, ha="left", va="top", fontsize=8, color="#777777")
            ax.set_xlabel("average soil relative humidity (%)")
            ax.set_ylabel(label)
            ax.set_title(label, loc="left")
            clean_axes(ax)
        fig.suptitle("Alpha diversity across the humidity gradient.", x=0.01, ha="left", fontsize=12)
        obs_rho, obs_p = stats.spearmanr(alpha_plot["average_soil_relative_humidity"], alpha_plot["observed_asvs"])
        shan_rho, shan_p = stats.spearmanr(alpha_plot["average_soil_relative_humidity"], alpha_plot["shannon_diversity"])
        if abs(shan_rho) > abs(obs_rho):
            caption = f"Shannon diversity shows the clearer humidity pattern here (rho={shan_rho:.2f}, p={fmt_p(shan_p)}); observed ASV richness is noisier (rho={obs_rho:.2f}, p={fmt_p(obs_p)})."
        else:
            caption = f"Observed ASV richness shows the clearer humidity pattern here (rho={obs_rho:.2f}, p={fmt_p(obs_p)}); Shannon diversity is rho={shan_rho:.2f}, p={fmt_p(shan_p)}."
        add_caption(fig, caption)
        plt.tight_layout()
        plt.show()
        '''
    )


def code_cell_association_plot() -> str:
    return dedent(
        '''
        count_cols = [col for col in counts_retained.columns if col.startswith("Atacama_ASV_")]
        q_counts = counts_retained[["sample_id", *count_cols]].merge(metadata[["sample_id", "average_soil_relative_humidity", "vegetation"]], on="sample_id")
        raw_counts = q_counts[count_cols].astype(float)
        # CLR uses a 0.5 pseudo-count so zero counts can be logged without creating infinite values.
        logged = np.log(raw_counts + 0.5)
        clr = logged.sub(logged.mean(axis=1), axis=0)
        sample_totals = counts_retained.set_index("sample_id")["total_reads"].reindex(q_counts["sample_id"]).to_numpy(dtype=float)
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
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.5), sharey=True)
        rng = np.random.default_rng(7)
        for ax, col, title in zip(axes, ["Raw p-value", "BH q-value"], ["Raw p-values", "BH q-values"]):
            y = humidity_results[col].to_numpy(dtype=float)
            x = rng.normal(0, 0.035, size=len(y))
            ax.scatter(x, y, s=18, color=OKABE_ITO["blue"], alpha=0.62, edgecolor="none")
            ax.axhline(0.05, color="#666666", ls="--", lw=1)
            ax.text(0.17, 0.05, "p or q = 0.05", ha="right", va="bottom", fontsize=8, color="#666666")
            ax.set_xlim(-0.18, 0.18)
            ax.set_xticks([])
            ax.set_title(title, loc="left", fontsize=10)
            clean_axes(ax)
        axes[0].set_ylabel("value")
        fig.suptitle("Raw p-values compared with BH q-values.", x=0.02, ha="left", fontsize=11)
        add_caption(fig, "Same 37 ASVs, before (left) and after (right) BH correction. Multiple-testing correction is conservative on purpose — it is the cost of asking many questions at once.")
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
                    arrow = "↑ with humidity" if row["Effect"] > 0 else "↓ with humidity"
                    effect_label = f"{row['Effect']:.2f} ({arrow})"
                else:
                    arrow = "↑ with vegetation" if row["Effect"] > 0 else "↓ with vegetation"
                    effect_label = f"{row['Effect']:.2f} ({arrow})"
                table_rows.append({
                    "Variable": variable,
                    "ASV": row["ASV"],
                    "Closest taxonomic match": row["closest_taxonomic_match"],
                    "Mean abundance (%)": round(row["Mean abundance (%)"], 2),
                    "Effect (Spearman rho or log2FC)": effect_label,
                    "Raw p-value": fmt_p(row["Raw p-value"]),
                    "BH q-value": fmt_p(row["BH q-value"]),
                    "Significant?": "✓" if row["BH q-value"] < 0.05 else "—",
                })

        q_table = pd.DataFrame(table_rows)

        def significant_style(row):
            if row["Significant?"] == "✓":
                return ["border-left: 3px solid #009E73"] + [""] * (len(row) - 1)
            return [""] * len(row)

        def sig_color(value):
            if value == "✓":
                return "color: #009E73; font-weight: 700"
            if value == "—":
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
                {"selector": "th", "props": [("text-align", "left"), ("background-color", "#f2f2f2"), ("border-bottom", "1px solid #999"), ("padding", "6px 7px")]},
                {"selector": "td", "props": [("padding", "6px 7px"), ("border-bottom", "1px solid #e6e6e6")]},
                {"selector": "tbody tr:nth-child(odd)", "props": [("background-color", "#fafafa")]},
            ])
        )
        '''
    )


def code_cell_lollipop() -> str:
    return dedent(
        '''
        significant = association_results.query("`BH q-value` < 0.05").copy()
        panels = ["Humidity", "Vegetation"]
        max_rows = int(significant.groupby("Variable").size().max()) if not significant.empty else 1
        fig, axes = plt.subplots(1, 2, figsize=(11.5, max(3.6, 0.42 * max(1, max_rows) + 1.5)))
        for ax, variable in zip(axes, panels):
            subset = significant.query("Variable == @variable").copy()
            if subset.empty:
                target = "humidity" if variable == "Humidity" else "vegetation"
                ax.text(0.5, 0.55, f"No ASVs reached q < 0.05\\nfor {target} in the QC'd dataset.", ha="center", va="center", transform=ax.transAxes, color="#666666", fontsize=10)
                ax.axvline(0, color="#777777", lw=1)
                ax.set_yticks([])
                ax.set_xlabel("Spearman rho" if variable == "Humidity" else "log2 fold-change")
                ax.set_title(variable, loc="left")
                clean_axes(ax)
                continue
            subset = subset.sort_values("Effect")
            y = np.arange(len(subset))
            colors = [phylum_color.get(row["phylum"], OKABE_ITO["blue"]) for _, row in subset.iterrows()]
            sizes = 45 + 18 * np.sqrt(subset["Mean abundance (%)"].clip(lower=0.01))
            ax.axvline(0, color="#777777", lw=1)
            ax.hlines(y, 0, subset["Effect"], color="#777777", lw=1.1)
            ax.scatter(subset["Effect"], y, s=sizes, color=colors, alpha=0.88, edgecolor="white", linewidth=0.4, zorder=3)
            labels = [asv_tax_label(row["ASV"]) for _, row in subset.iterrows()]
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=8)
            x_right = max(subset["Effect"].max(), 0) + 0.08
            for yi, (_, row) in enumerate(subset.iterrows()):
                ax.text(x_right, yi, f"q={fmt_p(row['BH q-value'])}", ha="left", va="center", fontsize=8, color="#999999")
            ax.set_xlabel("Spearman rho" if variable == "Humidity" else "log2 fold-change")
            ax.set_title(variable, loc="left")
            clean_axes(ax)
            ax.margins(x=0.22)
        fig.suptitle("ASVs with BH q-values below 0.05.", x=0.01, y=0.98, ha="left", fontsize=12)
        add_caption(fig, "Each bar is one ASV. Bars to the right are more abundant in higher-humidity (or vegetated) samples; bars to the left, lower-humidity (or non-vegetated). Dot size = overall abundance. Only ASVs with BH q < 0.05 are shown.")
        plt.tight_layout(rect=[0, 0.06, 1, 0.91])
        plt.show()
        '''
    )


def build_notebook(cache_paths: CachePaths) -> None:
    cache_files = {
        cache_paths.metadata.name: cache_paths.metadata.read_text(encoding="utf-8"),
        cache_paths.counts_retained.name: cache_paths.counts_retained.read_text(encoding="utf-8"),
        cache_paths.relative_top20.name: cache_paths.relative_top20.read_text(encoding="utf-8"),
        cache_paths.alpha.name: cache_paths.alpha.read_text(encoding="utf-8"),
        cache_paths.feature_key.name: cache_paths.feature_key.read_text(encoding="utf-8"),
        cache_paths.sequences_retained.name: cache_paths.sequences_retained.read_text(encoding="utf-8"),
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
    darwin_intro_image_uri = image_data_uri(DARWIN_INTRO_IMAGE_PATH)

    cells.append(
        md(
            dedent(
                f"""
                # Reading phylogenetic trees

                ## Section 1: What is a phylogenetic tree?

                A phylogenetic tree is a scientific hypothesis about relatedness. It is a branching diagram of descent: who shares ancestry with whom.

                The point of a tree is not to rank organisms from "simple" to "advanced." The point is to ask a sharper question: **which lineages share a more recent common ancestor?**

                Darwin sketched one of the most famous early tree diagrams in 1837. The small words above it were "I think" - a useful reminder that a tree is an evidence-based idea, not a decorative picture.

                <div style="display:flex; gap:18px; align-items:flex-start; margin:8px 0 10px 0;">
                  <div>
                    <img src="{darwin_intro_image_uri}" alt="Darwin I think tree sketch" width="390" style="max-width:100%; border:1px solid #e2e2e2;">
                    <div style="font-size:11px; color:#666; line-height:1.3; max-width:390px;">Darwin's Notebook B tree sketch, 1837. Public domain image via Wikimedia Commons; manuscript held by Cambridge University Library. Embedded here from the local notebook cache so it displays without internet image loading.</div>
                  </div>
                  <div style="max-width:520px; font-size:14px; line-height:1.5;">
                    <b>Tree-reading rule for today:</b><br>
                    More closely related means sharing a more recent common ancestor.<br><br>
                    This first section teaches three habits from tree-thinking research:
                    <ol>
                      <li>Read relationships by tracing branches to nodes.</li>
                      <li>Do not read across the tips like a left-to-right list.</li>
                      <li>Check what branch length means before interpreting distance.</li>
                    </ol>
                    The parts of a tree are:
                    <ul>
                      <li><b>Tip</b>: an organism or sequence being compared.</li>
                      <li><b>Branch</b>: one line of descent.</li>
                      <li><b>Internal node</b>: an inferred common ancestor.</li>
                      <li><b>Sister taxa</b>: two tips or groups that share an immediate common ancestor.</li>
                      <li><b>MRCA</b>: most recent common ancestor.</li>
                    </ul>
                  </div>
                </div>

                Tree-thinking approach adapted from Baum, Smith, and Donovan (2005), *The Tree-Thinking Challenge*, Science 310:979-980. The later DNA-to-tree flow is also guided by HHMI BioInteractive's phylogenetic-tree teaching materials.
                """
            )
        )
    )
    cells.append(code(make_setup_cell(cache_files) + "\n" + make_tree_helpers_cell() + "\n" + code_cell_annotated_mammal_tree()))
    cells.append(
        md(
            dedent(
                """
                ## Section 2: How to read relatedness

                To compare two tips, trace each one backward along the branches. Where their paths meet is their most recent common ancestor.

                If two tips share a more recent common ancestor with each other than with any other tip, they are closer relatives on that tree.

                Think: In the mammal tree, are dog and wolf closer to each other, or is dog closer to fox?
                """
            )
        )
    )
    cells.append(code(code_cell_tree_layouts()))
    cells.append(
        md(
            dedent(
                """
                ## Section 3: The tip-order trap

                The answer stayed the same in all three layouts: dog and wolf are sister taxa. The tree can look different because drawings can be rearranged, but the branching pattern still carries the relationship.

                A common mistake is to read across the page and assume nearby names are closest relatives. That is not enough. The better question is: **which tips share the most recent common ancestor?**
                """
            )
        )
    )
    cells.append(code(code_cell_rotation_diagram()))
    cells.append(
        md(
            dedent(
                """
                ## Section 4: Node rotation and branch length

                A tree can be rotated around an internal node without changing the relationships. Rotation changes the drawing, not the ancestry hypothesis.

                Branch length must be read from the tree's scale or caption. In the microbial tree later in this notebook, branch length means **DNA sequence difference in the 16S V4 region**, not time.

                Careful tree-reading now gives us the skill we need for real DNA data: tips, branches, internal nodes, sister groups, and MRCA.
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 5: From trees to real soil DNA

                Now that we know how to read a tree, we can use that skill on real biological data.

                Scientists collected soil samples from the Atacama Desert in Chile, one of the driest places on Earth. Even there, microbial life persists in the soil.

                The rest of this notebook asks: **which microbial DNA sequence types are present, what known groups do they most closely resemble, where are they abundant, and how are they related?**
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 6: Atacama dataset story

                This notebook uses **real QIIME 2 Atacama soil 16S data** from the 10% tutorial subsample: `data.qiime2.org/2024.10/tutorials/atacama-soils/10p/`. Soil was sampled across a humidity and aridity gradient from 22 sites along the Yungay and Baquedano transects, with metadata recording whether vegetation was present near each sample.

                Twelve samples with fewer than 100 reads after DADA2 denoising were excluded - at that depth, individual sequence counts are unreliable. This is standard quality control in microbiome studies. 46 samples remain across the humidity gradient.

                After QC, our humidity gradient is uneven - most surviving samples are from wetter sites, because drier sites yielded too little DNA to sequence reliably. This is the most common kind of sampling bias in extreme-environment microbiology.

                <div style="border:1px solid #d8d8d8; background:#f7f7f7; padding:10px 12px; width:600px; font-size:12px; line-height:1.35; color:#333;">
                Data source: real QIIME 2 Atacama Soils 10p tutorial subsample (data.qiime2.org/2024.10/tutorials/atacama-soils/). Processed with QIIME 2 Amplicon 2024.10.1, DADA2 denoising, SILVA 138 99% Naive Bayes classifier (sha256 c08a1aa4...62613616). 401 ASVs across 61 samples in the raw output; 37 ASVs across 46 samples after QC.
                </div>
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 7: What is a 16S ASV?

                Scientists often study bacteria by sequencing a marker gene called **16S rRNA**. After the reads are cleaned by DADA2, the workflow identifies exact DNA sequence patterns called **ASVs**, or Amplicon Sequence Variants.

                An ASV is not automatically a species. It is a precise 16S sequence pattern found in the samples.

                Keep three ideas separate: **abundance** means how much of an ASV is in a sample; **taxonomic match** means what known group the sequence resembles in SILVA; **sequence relatedness** means how ASV sequences cluster in a tree.

                Many real environmental ASVs do not have known genus or species labels. You will see "Unassigned" or family-level labels on some tips, and that is normal: much microbial diversity is still unnamed.
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 8: Load data and apply QC

                The notebook starts from real Atacama tables embedded with the notebook so every student sees the same data when pressing Run all. We apply the same QC rules before any tree or statistics: keep samples with enough reads, then keep ASVs seen in at least three usable samples.
                """
            )
        )
    )
    cells.append(code(code_cell_dataset_summary()))
    cells.append(
        md(
            dedent(
                """
                *The dataset is sparse, and that is part of the lesson: QC removes unreliable samples before we make any biological claim.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 9: Alignment of representative ASV sequences

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
                *ASVs that look similar across the bright columns are more closely related. The phylogenetic tree (section 11) formalizes this intuition.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 10: Distance matrix

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
                ## Section 11: UPGMA tree (the payoff)

                Tips are ASVs. Branch length shows sequence difference in this 252-bp V4 region. Branches that meet recently, with a short path between them, suggest closer sequence relatedness - our best evidence of evolutionary relatedness, but not proof of it.
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
                ## Section 12: Relative abundance - what is actually in these samples?

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
                ## Section 13: Alpha diversity

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
                *Use the trend lines as a visual guide, then describe the direction carefully. Real shallow data can show modest correlations rather than dramatic ones.*
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Section 14: BH-corrected association tests

                We test 37 ASVs to see if abundance changes with humidity. With a p < 0.05 cutoff, we would expect about 37 × 0.05 ≈ 2 false positives by random chance alone, even if nothing is truly associated. Benjamini-Hochberg (BH) correction adjusts for this. A q-value of 0.05 means we expect about 5% of discoveries below that threshold to be false alarms.

                These q-values belong to abundance-versus-metadata tests. They are not tree branch support values, and they do not prove that an ASV is a species.
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
                ## Section 15: Final student report

                Fill in the report using the section numbers named in each question.

                1. Which 3 ASVs are most abundant in your samples? Refer to section 12 table.
                   [your answer]

                2. Which ASVs are significantly associated with humidity? Refer to section 14 lollipop.
                   [your answer]

                3. Which ASVs are significantly associated with vegetation? Refer to section 14 lollipop.
                   [your answer]

                4. What does alpha diversity suggest about humid versus arid samples? Refer to section 13.
                   [your answer]

                5. In the UPGMA tree, which ASVs cluster closest together? Refer to section 11.
                   [your answer]

                6. Are closely related ASVs from section 11 also similar in abundance from section 12 or in humidity association from section 14?
                   [your answer]

                7. Synthesis: An ASV can be (a) very abundant, (b) statistically associated with humidity, and (c) closely related to another ASV in the tree - and these are three different things. Explain in 2-3 sentences why these are three different ideas and why all three matter for understanding the soil microbiome.
                   [your answer]
                """
            )
        )
    )
    cells.append(
        md(
            dedent(
                """
                ## Appendix: Full data provenance

                - Source URL: `data.qiime2.org/2024.10/tutorials/atacama-soils/10p/`
                - QIIME 2 version: Amplicon 2024.10.1, run in WSL
                - DADA2 denoising parameters: defaults from the QIIME 2 Atacama tutorial
                - SILVA classifier: SILVA 138 99% OTUs full-length Naive Bayes classifier, sha256 `c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616`
                - Raw artifact counts: 401 ASVs, 61 samples
                - After sample QC (`>=100` reads plus complete humidity and vegetation metadata): 46 samples
                - After ASV prevalence filter (`>=3` samples): 37 ASVs
                - Study reference: Neilson et al. 2017, mSystems, https://doi.org/10.1128/mSystems.00195-16
                """
            )
        )
    )

    nb["cells"] = cells
    if not (28 <= len(cells) <= 36):
        raise AssertionError(f"Goal 2 requires 28-36 cells; built {len(cells)}")
    NOTEBOOK_PATH.write_text(nbf.writes(nb), encoding="utf-8", newline="\n")


def main() -> None:
    cache_paths = build_cache()
    build_notebook(cache_paths)
    print(f"Wrote {NOTEBOOK_PATH.name}")
    print(f"Wrote Goal 2 cache files under {CACHE_DIR}")


if __name__ == "__main__":
    main()
