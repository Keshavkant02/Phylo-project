from __future__ import annotations

import csv
import json
import math
import tempfile
import zipfile
from collections import OrderedDict
from pathlib import Path

import h5py
import numpy as np
from scipy import stats


ROOT = Path(__file__).resolve().parent
SOURCE_DIR = ROOT / "tmp" / "atacama_qiime2_source"
CACHE_DIR = ROOT / "soil_16s_class_cache"
TOP_N = 12

METADATA_PATH = SOURCE_DIR / "sample_metadata.tsv"
TABLE_QZA = SOURCE_DIR / "atacama-table.qza"
REP_SEQS_QZA = SOURCE_DIR / "atacama-rep-seqs.qza"


def read_tsv_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            sample_id = row["sample-id"]
            if sample_id == "#q2:types":
                continue
            rows[sample_id] = row
    return rows


def read_qza_file(qza_path: Path, suffix: str) -> bytes:
    with zipfile.ZipFile(qza_path) as zf:
        matches = [name for name in zf.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f"Expected one {suffix} in {qza_path}, found {matches}")
        return zf.read(matches[0])


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
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = []
        elif current is not None:
            records[current].append(line.strip())
    return {key: "".join(chunks) for key, chunks in records.items()}


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def bh_adjust(p_values: list[float]) -> list[float]:
    n = len(p_values)
    order = sorted(range(n), key=lambda idx: p_values[idx])
    adjusted = [math.nan] * n
    running = 1.0
    for rank_from_end, idx in enumerate(reversed(order), start=1):
        rank = n - rank_from_end + 1
        value = min(running, p_values[idx] * n / rank)
        running = value
        adjusted[idx] = min(value, 1.0)
    return adjusted


def shannon(counts: np.ndarray) -> float:
    total = float(counts.sum())
    if total <= 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(-(p * np.log(p)).sum())


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def main() -> None:
    for path in [METADATA_PATH, TABLE_QZA, REP_SEQS_QZA]:
        if not path.exists() or path.stat().st_size < 100:
            raise FileNotFoundError(f"Missing or invalid source artifact: {path}")

    CACHE_DIR.mkdir(exist_ok=True)

    metadata = read_tsv_metadata(METADATA_PATH)
    feature_ids, sample_ids, table = read_biom_table(TABLE_QZA)
    sequences = read_fasta_from_qza(REP_SEQS_QZA)

    missing_metadata = sorted(set(sample_ids) - set(metadata))
    if missing_metadata:
        raise ValueError(f"Feature table samples missing metadata: {missing_metadata[:5]}")

    total_by_feature = table.sum(axis=1)
    prevalence_by_feature = (table > 0).sum(axis=1)
    top_indices = list(np.argsort(-total_by_feature)[:TOP_N])
    top_features = [feature_ids[idx] for idx in top_indices]
    asv_labels = {feature_id: f"Atacama_ASV_{i:02d}" for i, feature_id in enumerate(top_features, start=1)}

    sample_rows: list[dict[str, object]] = []
    humidity = []
    vegetation = []
    for sample_id in sample_ids:
        row = metadata[sample_id]
        humidity_value = safe_float(row["average-soil-relative-humidity"])
        humidity.append(humidity_value)
        vegetation.append(row["vegetation"])
        sample_rows.append(
            {
                "sample_id": sample_id,
                "transect_name": row["transect-name"],
                "site_name": row["site-name"],
                "depth": row["depth"],
                "elevation": row["elevation"],
                "average_soil_relative_humidity": row["average-soil-relative-humidity"],
                "vegetation": row["vegetation"],
                "ph": row["ph"],
                "toc": row["toc"],
                "ec": row["ec"],
                "percentcover": row["percentcover"],
            }
        )

    write_csv(
        CACHE_DIR / "atacama_sample_metadata_mini.csv",
        sample_rows,
        [
            "sample_id",
            "transect_name",
            "site_name",
            "depth",
            "elevation",
            "average_soil_relative_humidity",
            "vegetation",
            "ph",
            "toc",
            "ec",
            "percentcover",
        ],
    )

    top_table = table[top_indices, :].T
    sample_totals = table.sum(axis=0)
    top_count_rows = []
    rel_rows = []
    for sample_position, sample_id in enumerate(sample_ids):
        count_row = {"sample_id": sample_id}
        rel_row = {"sample_id": sample_id}
        top_sum = 0.0
        for feature_id, count in zip(top_features, top_table[sample_position]):
            label = asv_labels[feature_id]
            count_row[label] = int(count)
            rel_value = 100.0 * float(count) / sample_totals[sample_position] if sample_totals[sample_position] else 0.0
            rel_row[label] = round(rel_value, 5)
            top_sum += float(count)
        count_row["other_asvs"] = int(sample_totals[sample_position] - top_sum)
        count_row["total_reads"] = int(sample_totals[sample_position])
        rel_row["other_asvs"] = (
            round(100.0 * (sample_totals[sample_position] - top_sum) / sample_totals[sample_position], 5)
            if sample_totals[sample_position]
            else 0.0
        )
        rel_row["total_reads"] = int(sample_totals[sample_position])
        top_count_rows.append(count_row)
        rel_rows.append(rel_row)

    asv_columns = [asv_labels[feature_id] for feature_id in top_features]
    write_csv(CACHE_DIR / "atacama_feature_table_top12.csv", top_count_rows, ["sample_id", *asv_columns, "other_asvs", "total_reads"])
    write_csv(CACHE_DIR / "atacama_relative_abundance_top12.csv", rel_rows, ["sample_id", *asv_columns, "other_asvs", "total_reads"])

    feature_rows = []
    for idx, feature_id in zip(top_indices, top_features):
        feature_rows.append(
            {
                "asv_label": asv_labels[feature_id],
                "qiime_feature_id": feature_id,
                "total_reads": int(total_by_feature[idx]),
                "prevalence_samples": int(prevalence_by_feature[idx]),
                "representative_sequence_length": len(sequences.get(feature_id, "")),
            }
        )
    write_csv(
        CACHE_DIR / "atacama_feature_key.csv",
        feature_rows,
        ["asv_label", "qiime_feature_id", "total_reads", "prevalence_samples", "representative_sequence_length"],
    )

    fasta_lines = []
    for feature_id in top_features:
        label = asv_labels[feature_id]
        fasta_lines.append(f">{label} qiime_feature_id={feature_id}")
        sequence = sequences.get(feature_id, "")
        fasta_lines.extend(sequence[i : i + 80] for i in range(0, len(sequence), 80))
    (CACHE_DIR / "atacama_top_asv_sequences.fasta").write_text("\n".join(fasta_lines).strip() + "\n", encoding="utf-8", newline="\n")

    sample_has_reads = sample_totals > 0
    rel_matrix = np.zeros_like(top_table, dtype=float)
    rel_matrix[sample_has_reads, :] = top_table[sample_has_reads, :] / sample_totals[sample_has_reads, None]
    humidity_array = np.asarray(humidity, dtype=float)
    valid_humidity = np.isfinite(humidity_array) & sample_has_reads
    veg_yes = np.asarray([value == "yes" for value in vegetation])
    veg_no = np.asarray([value == "no" for value in vegetation])

    stat_rows = []
    spearman_p = []
    vegetation_p = []
    for col_index, feature_id in enumerate(top_features):
        rel = rel_matrix[:, col_index]
        rho, p_value = stats.spearmanr(humidity_array[valid_humidity], rel[valid_humidity])
        if np.isnan(rho) or np.isnan(p_value):
            rho, p_value = 0.0, 1.0
        try:
            mw = stats.mannwhitneyu(rel[veg_yes & sample_has_reads], rel[veg_no & sample_has_reads], alternative="two-sided")
            vegetation_p_value = float(mw.pvalue)
        except ValueError:
            vegetation_p_value = 1.0
        spearman_p.append(float(p_value))
        vegetation_p.append(vegetation_p_value)
        stat_rows.append(
            {
                "asv_label": asv_labels[feature_id],
                "spearman_rho_vs_humidity": round(float(rho), 5),
                "spearman_p_vs_humidity": float(p_value),
                "mannwhitney_p_vegetated_vs_unvegetated": vegetation_p_value,
                "mean_relative_abundance_percent": round(float(np.mean(rel[sample_has_reads]) * 100), 5),
                "max_relative_abundance_percent": round(float(np.max(rel[sample_has_reads]) * 100), 5),
            }
        )

    spearman_q = bh_adjust(spearman_p)
    vegetation_q = bh_adjust(vegetation_p)
    for row, q_humidity, q_vegetation in zip(stat_rows, spearman_q, vegetation_q):
        row["spearman_q_bh_vs_humidity"] = q_humidity
        row["mannwhitney_q_bh_vegetated_vs_unvegetated"] = q_vegetation

    write_csv(
        CACHE_DIR / "atacama_top_asv_stats.csv",
        stat_rows,
        [
            "asv_label",
            "mean_relative_abundance_percent",
            "max_relative_abundance_percent",
            "spearman_rho_vs_humidity",
            "spearman_p_vs_humidity",
            "spearman_q_bh_vs_humidity",
            "mannwhitney_p_vegetated_vs_unvegetated",
            "mannwhitney_q_bh_vegetated_vs_unvegetated",
        ],
    )

    alpha_rows = []
    observed = (table > 0).sum(axis=0)
    shannon_values = np.asarray([shannon(table[:, sample_index]) for sample_index in range(table.shape[1])])
    for sample_position, sample_id in enumerate(sample_ids):
        alpha_rows.append(
            {
                "sample_id": sample_id,
                "total_reads": int(sample_totals[sample_position]),
                "observed_asvs": int(observed[sample_position]),
                "shannon_entropy": round(float(shannon_values[sample_position]), 5),
            }
        )
    write_csv(CACHE_DIR / "atacama_alpha_diversity.csv", alpha_rows, ["sample_id", "total_reads", "observed_asvs", "shannon_entropy"])

    alpha_stat_rows = []
    alpha_p = []
    for metric, values in [
        ("total_reads", sample_totals),
        ("observed_asvs", observed.astype(float)),
        ("shannon_entropy", shannon_values),
    ]:
        rho, p_value = stats.spearmanr(humidity_array[valid_humidity], values[valid_humidity])
        if np.isnan(rho) or np.isnan(p_value):
            rho, p_value = 0.0, 1.0
        alpha_stat_rows.append(
            {
                "metric": metric,
                "spearman_rho_vs_humidity": round(float(rho), 5),
                "spearman_p_vs_humidity": float(p_value),
            }
        )
        alpha_p.append(float(p_value))
    for row, q_value in zip(alpha_stat_rows, bh_adjust(alpha_p)):
        row["spearman_q_bh_vs_humidity"] = q_value
    write_csv(
        CACHE_DIR / "atacama_alpha_diversity_stats.csv",
        alpha_stat_rows,
        ["metric", "spearman_rho_vs_humidity", "spearman_p_vs_humidity", "spearman_q_bh_vs_humidity"],
    )

    manifest = {
        "title": "Atacama soil microbiome mini-cache for browser-only teaching",
        "source_tutorial": "QIIME 2 2024.10 Atacama soil microbiome tutorial and q2-vsearch chimera tutorial",
        "source_study": "Neilson et al. 2017 mSystems, Significant Impacts of Increasing Aridity on the Arid Soil Microbiome",
        "source_doi": "10.1128/mSystems.00195-16",
        "source_urls": {
            "sample_metadata": "https://data.qiime2.org/2024.10/tutorials/atacama-soils/sample_metadata.tsv",
            "feature_table": "https://data.qiime2.org/2024.10/tutorials/chimera/atacama-table.qza",
            "representative_sequences": "https://data.qiime2.org/2024.10/tutorials/chimera/atacama-rep-seqs.qza",
        },
        "sample_count": len(sample_ids),
        "feature_count_full_table": len(feature_ids),
        "top_feature_count_in_cache": TOP_N,
        "teaching_use": "Use for real soil ASV abundance, humidity/vegetation association, and adjusted p-value teaching without installing QIIME 2.",
    }
    (CACHE_DIR / "atacama_mini_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")

    readme = """# Atacama Soil Mini-Cache

This folder includes a small derived cache from the QIIME 2 Atacama soil microbiome tutorial.

Use this for browser-only teaching of real soil ASV abundance patterns. The heavy QIIME 2 steps are not run in Colab; this cache stores the derived tables students need.

Files:

- `atacama_sample_metadata_mini.csv`: selected sample metadata for table samples.
- `atacama_feature_table_top12.csv`: counts for the 12 most abundant ASVs plus all other ASVs.
- `atacama_relative_abundance_top12.csv`: relative abundance percentages for the same ASVs.
- `atacama_feature_key.csv`: ASV label to original QIIME feature ID mapping.
- `atacama_top_asv_sequences.fasta`: representative sequences for top ASVs.
- `atacama_top_asv_stats.csv`: Spearman humidity and Mann-Whitney vegetation tests with Benjamini-Hochberg q-values.
- `atacama_alpha_diversity.csv`: total reads, observed ASVs, and Shannon entropy per sample.
- `atacama_alpha_diversity_stats.csv`: alpha-diversity correlation with humidity and BH q-values.
- `atacama_mini_manifest.json`: provenance and source URLs.

Teaching caveat: adjusted p-values here test abundance/metadata associations. They are not tree branch support values.
"""
    (CACHE_DIR / "ATACAMA_MINI_README.md").write_text(readme, encoding="utf-8", newline="\n")

    print(json.dumps({"status": "passed", "samples": len(sample_ids), "features": len(feature_ids), "top_features": TOP_N}, indent=2))


if __name__ == "__main__":
    main()
