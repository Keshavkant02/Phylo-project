from __future__ import annotations

import csv
import heapq
import io
import math
import re
import zipfile
from collections import Counter, OrderedDict
from dataclasses import dataclass
from pathlib import Path

from Bio.Align import PairwiseAligner


ROOT = Path(__file__).resolve().parent
CACHE_DIR = ROOT / "soil_16s_class_cache"
TMP_DIR = ROOT / "tmp"

QUERY_FASTA = CACHE_DIR / "goal2_atacama_rep_seqs_retained_asvs.fasta"
FEATURE_KEY = CACHE_DIR / "goal2_atacama_feature_key.csv"
OUT_ASSIGNMENTS = CACHE_DIR / "goal2_atacama_silva_static_taxonomy_assignments.csv"

SILVA_SEQS_QZA = TMP_DIR / "silva-138-99-seqs-515-806.qza"
SILVA_TAX_QZA = TMP_DIR / "silva-138-99-tax-515-806.qza"

KMER = 11
TOP_KMER_CANDIDATES = 350
TOP_REPORT_HITS = 8


def read_qza_text(qza_path: Path, suffix: str) -> str:
    if not qza_path.exists() or not zipfile.is_zipfile(qza_path):
        raise FileNotFoundError(f"Missing valid QIIME artifact: {qza_path}")
    with zipfile.ZipFile(qza_path) as zf:
        matches = [name for name in zf.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f"Expected one {suffix} in {qza_path}, found {matches}")
        return zf.read(matches[0]).decode("utf-8")


def parse_fasta(text: str) -> OrderedDict[str, str]:
    records: OrderedDict[str, list[str]] = OrderedDict()
    current = None
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = []
        elif current is not None:
            records[current].append(re.sub("[^ACGTN-]", "", line.strip().upper()))
    return OrderedDict((key, "".join(chunks).replace("-", "")) for key, chunks in records.items())


def read_fasta(path: Path) -> OrderedDict[str, str]:
    return parse_fasta(path.read_text(encoding="utf-8"))


def parse_taxonomy_qza(path: Path) -> dict[str, str]:
    text = read_qza_text(path, "taxonomy.tsv")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    out: dict[str, str] = {}
    for row in reader:
        feature_id = row.get("Feature ID") or row.get("feature-id") or row.get("id")
        taxon = row.get("Taxon") or row.get("taxonomy") or ""
        if feature_id and taxon:
            out[feature_id] = taxon
    return out


def clean_rank_name(value: str) -> str:
    value = value.strip()
    value = re.sub(r"^[dkpcofgs]__", "", value)
    value = re.sub(r"^D_[0-6]__", "", value)
    if not value or value.lower() in {
        "uncultured",
        "unidentified",
        "unknown",
        "metagenome",
        "unassigned",
        "uncultured bacterium",
    }:
        return ""
    return value


def parse_ranks(taxon: str) -> dict[str, str]:
    chunks = [chunk.strip() for chunk in taxon.split(";")]
    ranks = {"phylum": "", "family": "", "genus": ""}
    for chunk in chunks:
        if chunk.startswith(("p__", "D_1__")):
            ranks["phylum"] = clean_rank_name(chunk)
        elif chunk.startswith(("f__", "D_4__")):
            ranks["family"] = clean_rank_name(chunk)
        elif chunk.startswith(("g__", "D_5__")):
            genus = clean_rank_name(chunk)
            if genus.endswith("aceae") or genus.endswith("ales"):
                # SILVA sometimes repeats a higher-rank placeholder at genus level.
                # Keep these as family/order-level labels instead of presenting them as genera.
                genus = ""
            ranks["genus"] = genus
    return ranks


def display_match(ranks: dict[str, str]) -> str:
    return ranks["genus"] or ranks["family"] or ranks["phylum"] or "Unassigned at genus level"


def kmers(seq: str, k: int = KMER) -> set[str]:
    seq = seq.upper().replace("N", "")
    if len(seq) < k:
        return set()
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


@dataclass(frozen=True)
class Candidate:
    ref_id: str
    shared_kmers: int


@dataclass(frozen=True)
class Hit:
    ref_id: str
    score: float
    identity: float
    coverage: float
    matches: int
    aligned: int
    shared_kmers: int
    taxon: str


def local_alignment_hit(query: str, ref: str, ref_id: str, shared_kmers: int, taxonomy: dict[str, str], aligner: PairwiseAligner) -> Hit:
    aln = aligner.align(query, ref)[0]
    matches = 0
    aligned = 0
    for (qs, qe), (rs, re) in zip(aln.aligned[0], aln.aligned[1]):
        q_block = query[qs:qe]
        r_block = ref[rs:re]
        for q_base, r_base in zip(q_block, r_block):
            if q_base in "ACGT" and r_base in "ACGT":
                aligned += 1
                matches += int(q_base == r_base)
    identity = matches / aligned if aligned else 0.0
    coverage = aligned / len(query) if query else 0.0
    return Hit(
        ref_id=ref_id,
        score=float(aln.score),
        identity=identity,
        coverage=coverage,
        matches=matches,
        aligned=aligned,
        shared_kmers=shared_kmers,
        taxon=taxonomy.get(ref_id, ""),
    )


def best_hits_for_queries(queries: OrderedDict[str, str], refs: OrderedDict[str, str], taxonomy: dict[str, str]) -> dict[str, list[Hit]]:
    query_kmers = {qid: kmers(seq) for qid, seq in queries.items()}
    kmer_to_queries: dict[str, list[str]] = {}
    for qid, words in query_kmers.items():
        for word in words:
            kmer_to_queries.setdefault(word, []).append(qid)

    heaps: dict[str, list[tuple[int, str]]] = {qid: [] for qid in queries}

    for ref_id, ref_seq in refs.items():
        counts: Counter[str] = Counter()
        for word in kmers(ref_seq):
            for qid in kmer_to_queries.get(word, []):
                counts[qid] += 1
        for qid, count in counts.items():
            heap = heaps[qid]
            item = (int(count), ref_id)
            if len(heap) < TOP_KMER_CANDIDATES:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5

    results: dict[str, list[Hit]] = {}
    for qid, query in queries.items():
        candidates = sorted((Candidate(ref_id=ref_id, shared_kmers=count) for count, ref_id in heaps[qid]), key=lambda x: (-x.shared_kmers, x.ref_id))
        refined = [
            local_alignment_hit(query, refs[candidate.ref_id], candidate.ref_id, candidate.shared_kmers, taxonomy, aligner)
            for candidate in candidates
        ]
        refined.sort(key=lambda h: (h.coverage >= 0.8, h.identity * h.coverage, h.identity, h.coverage, h.score), reverse=True)
        results[qid] = refined[:TOP_REPORT_HITS]
    return results


def load_feature_id_by_asv() -> dict[str, str]:
    with FEATURE_KEY.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return {row["asv"]: row["qiime_feature_id"] for row in reader}


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.4f}"


def write_assignments(hits_by_query: dict[str, list[Hit]]) -> None:
    feature_ids = load_feature_id_by_asv()
    fieldnames = [
        "asv",
        "qiime_feature_id",
        "closest_taxonomic_match",
        "phylum",
        "family",
        "genus",
        "top_silva_feature_id",
        "percent_identity",
        "query_coverage",
        "shared_11mers",
        "top_taxon",
        "method",
        "reference_source",
    ]
    rows = []
    for asv, hits in hits_by_query.items():
        hit = hits[0]
        ranks = parse_ranks(hit.taxon)
        rows.append(
            {
                "asv": asv,
                "qiime_feature_id": feature_ids.get(asv, ""),
                "closest_taxonomic_match": display_match(ranks),
                "phylum": ranks["phylum"] or "Unassigned",
                "family": ranks["family"],
                "genus": ranks["genus"],
                "top_silva_feature_id": hit.ref_id,
                "percent_identity": fmt(hit.identity * 100),
                "query_coverage": fmt(hit.coverage * 100),
                "shared_11mers": hit.shared_kmers,
                "top_taxon": hit.taxon,
                "method": "nearest_silva_515_806_reference_local_alignment_not_qiime_naive_bayes",
                "reference_source": "QIIME 2 2024.10 SILVA 138 SSURef NR99 515F/806R reference sequences and taxonomy",
            }
        )

    with OUT_ASSIGNMENTS.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    queries = read_fasta(QUERY_FASTA)
    refs = parse_fasta(read_qza_text(SILVA_SEQS_QZA, "dna-sequences.fasta"))
    taxonomy = parse_taxonomy_qza(SILVA_TAX_QZA)
    missing_tax = len(set(refs) - set(taxonomy))
    if missing_tax:
        raise ValueError(f"{missing_tax} SILVA references lack taxonomy labels")
    hits = best_hits_for_queries(queries, refs, taxonomy)
    write_assignments(hits)
    print(f"Wrote {OUT_ASSIGNMENTS.relative_to(ROOT)} with {len(hits)} ASV assignments against {len(refs)} SILVA references.")


if __name__ == "__main__":
    main()
