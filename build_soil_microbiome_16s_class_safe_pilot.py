from __future__ import annotations

import ast
import csv
import json
import math
import textwrap
from pathlib import Path
from xml.sax.saxutils import escape as xml_escape


ROOT = Path(__file__).resolve().parent
CACHE_DIR = ROOT / "soil_16s_class_cache"
NOTEBOOK_PATH = ROOT / "soil_microbiome_16s_class_safe_colab.ipynb"
RETRIEVED_DATE = "2026-05-25"


RAW_REFERENCE_FASTA = r""">NR_102783.2 Bacillus subtilis subsp. subtilis strain 168 16S ribosomal RNA, complete sequence
TTATCGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGGA
CAGATGGGAGCTTGCTCCCTGATGTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTGTAAGA
CTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATGGTTGTTTGAACCGCATGGTTCAAACATAAAAG
GTGGCTTCGGCTACCACTTACAGATGGACCCGCGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAA
GGCGACGATGCGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTA
CGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGA
AGGTTTTCGGATCGTAAAGCTCTGTTGTTAGGGAAGAACAAGTGCCGTTCGAATAGGGCGGTACCTTGAC
GGTACCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTG
TCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAAC
CGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTG
AAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGAG
CGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGT
TAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAG
ACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGC
GAAGAACCTTACCAGGTCTTGACATCCTCTGACAATCCTAGAGATAGGACGTCCCCTTCGGGGGCAGAGT
GACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAAC
CCTTGATCTTAGTTGCCAGCATTCAGTTGGGCACTCTAAGGTGACTGCCGGTGACAAACCGGAGGAAGGT
GGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACAGAACAAAG
GGCAGCGAAACCGCGAGGTTAAGCCAATCCCACAAATCTGTTCTCAGTTCGGATCGCAGTCTGCAACTCG
ACTGCGTGAAGCTGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGT
ACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCTTTTAGGAGCCAGC
CGCCGAAGGTGGGACAGATGATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGAT
CACCTCCTTT

>NR_115715.1 Pseudomonas fluorescens strain CCM 2115 16S ribosomal RNA, partial sequence
AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGTAGAGAGA
AGCTTGCTTCTCTTGAGAGCGGCGGACGGGTGAGTAAAGCCTAGGAATCTGCCTGGTAGTGGGGGATAAC
GTTCGGAAACGGACGCTAATACCGCATACGTCCTACGGGAGAAAGCAGGGGACCTTCGGGCCTTGCGCTA
TCAGATGAGCCTAGGTCGGATTAGCTAGTTGGTGAGGTAATGGCTCACCAAGGCGACGATCCGTAACTGG
TCTGAGAGGATGATCAGTCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGA
ATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAG
CACTTTAAGTTGGGAGGAAGGGCATTAACCTAATACGTTAGTGTTTTGACGTTACCGACAGAATAAGCAC
CGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAA
AGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAAC
TGACTGACTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAG
GAACACCAGTGGCGAAGGCGACCACCTGGACTAATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAAC
AGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCAACTAGCCGTTGGGAGCCTTGAGCTCTTAG
TGGCGCAGCTAACGCATTAAGTTGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGA
CGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGCCTTG
ACATCCAATGAACTTTCTAGAGATAGATTGGTGCCTTCGGGAACATTGAGACAGGTGCTGCATGGCTGTC
GTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGTAACGAGCGCAACCCTTGTCCTTAGTTACCAGCA
CGTAATGGTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCAT
CATGGCCCTTACGGCCTGGGCTACACACGTGCTACAATGGTCGGTACAGAGGGTTGCCAAGCCGCGAGGT
GGAGCTAATCCCACAAAACCGATCGTAGTCCGGATCGCAGTCTGCAACTCGACTGCGTGAAGTCGGAATC
GCTAGTAATCGCGAATCAGAATGTCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACC
ATGGGAGTGGGTTGCACCAGAAGTAGCTAGTCTAACCTTCGGGAGGACGGTTACCACGGTGTGATTCATG
ACTGGGGTGAAGTCGTAACAAGGTAGCCGTAGGGGAACCTGCGGCTGGAT

>Y00411.1 Streptomyces coelicolor 16S rRNA gene from rrnD
TGGGCCCGCATCACCATCGGCGTCCTCGCCGAGCTGGCCTTCCTGGCCTACGTCTACGTTCTGGGCGGCC
GAGCCGTGCGCGACGGCGAGACGGGTGACGTCGAGGCAGCCGAACGCAGCGCCACGGTGCCAACAGCCGC
CTGATGTGCATCCACCCCTGCGAGCTGCTAGTGTCCTCTTCGTTCCCGCAAGAGCCGTTGACACGGAGCG
AGCGGGGAGGTAGATTCGAACAGTTGCCTGGAGACGGGTTCACCCCAGAGGGCAACAGTGAACATCTACC
AGCTTCTCCGAATCAACGAATTCGACGAAGCACTCTCCCGATGAATCGGAAACGAAGGCCGGTAAGACCG
GCTCGAAAGTTCTGATAAAGTCGGAGCCGCCGGAAAGGGAAACGCGAAAGCGGGAACCTGGAAAGCGCCG
AGGAAATCGGATCGGAAAGATCTGATAGAGTCGGAAACGCAAGACCGAAGGGAAGCGCCCGGAGGAAAGC
CCGAGAGGGTGAGTACAAAGGAAGCGTCCGTTCCTTGAGAACTCAACAGCGTGCCAAAAGTCAACGCCAG
ATATGTTGATACCCCGACCTGATCGGATCTCCGTTCGGGTTGAGGTTCCTTTGAAGTAACACAACAGCGA
GGACGCTGTGAACGGTCGGATTATTCCTCCGACTGTTCCGCTCTCGTGGTGTCACCCGATTACGGGTATA
CATTCACGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACG
ATGAACCACTTCGGTGGGGATTAGTGGCGAACGGGTGAGTAACACGTGGGCAATCTGCCCTTCACTCTGG
GACAAGCCCTGGAAACGGGGTCTAATACCGGATACTGACCCTCGCAGGCATCTGCGAGGTTCGAAAGCTC
CGGCGGTGAAGGATGAGCCCGCGGCCTATCAGCTTGTTGGTGAGGTAATGGCTCACCAAGGCGACGACGG
GTAGCCGGCCTGAGAGGGCGACCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGC
AGTGGGGAATGTTGCACAATGGGCGAAAGCCTGATGCAGCGACGCCGCGTGAGGGATGACGGCCTTCGGG
TTGTAAACCTCTTTCAGCAGGGAAGAAGCGAAAGTGACGGTACCTGCAGAAGAAGCGCCGGCTAACTACG
TGCCAGCAGCCGCGGTAATACGTAGGGCGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGG
CGGCTTGTCACGTCGGTTGTGAAAGCCCGGGGCTTAACCCCGCCACTGCAGTCGATACGGGCAGGCTAGA
GTTCGGTAGGGGAGATCGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGGTGG
CGAAGGCGGATCTCTGGGCCGATACTGACGCTGAGGAGCGAAAGNGTGGGGAGCGAACAGGATTAGATAC
CCTGGTAGTCCACGCCGTAAACGGTGGGCACTAGGTGTGGGCAACATTCCACGTTGTCCGTGCCGCAGCT
AACGCATTAAGTGCCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGGAATTGACGGGGGCCCG
CACAAGCGGCGGAGCATGTGGCTTAATTCGACGCAACGCGAAGAACCTTACCAAGGCTTGACATACACCG
GAAAGCATCAGAGATGGTGCCCCCCTTGTGGTCGGTGTACAGGTGGTGCATGGCTGTCGTCAGCTCGTGT
CGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCCCGTGTTGCCAGCAAGCCCTTCGGGG
TGTTGGGGACTCACGGGAGACCGCCGGGGTCAACTCGGAGGAAGGTGGGGACGACGTCAAGTCATCATGC
CCCTTATGTCTTGGGCTGCACACGTGCTACAATGGCCGGTACAATGAGCTGCGATACCGCAAGGTGGAGC
GAATCTCAAAAAGCCGGTCTCAGTTCGGATTGGGGTCTGCAACTCGACCCCATGAAGTCGGAGTCGCTAG
TAATCGCAGATCAGCATTGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGTCACGA
AAGTCGGTAACACCCGAAGCCGGTGGCCCAACCCCTTGTGGGAGGGAGCTGTCGAAGGTGGGACTGGCGA
TTGGGACGAAGTCGTAACAAGGTAGCCGTACCGGAAGGTGCGGCTGGATCACCTCCTTTCTAAGGAGCAC
ATAGCCGACTGCAGCGAAATGTCCTGCACGGTTGCTCATGGGTGGAACGTTGACTACTCGGCACGGTCTT
CTTGATGGATCACTAGTACTGCTTCGGCGTGGAACGTGACTTCAAAGAGGGGTTCGTGTCGGGCACGCTG
TTGGGTATCTGAGGGTACGGCCGTGAGGTCGCCTTCAGTTGCCGGCCCCGGTAAAAATCCGCGTGAGTGG
GTTGTGACGGGTGGTTGGTCGTTGTTTGAGAACTGCACAGTGGACGCGAGCATCTGTGGCCAAGTTTTTA
AGGGCGCACGGTGGATGCCTT

>D14513.1 Rhizobium leguminosarum gene for 16S rRNA, complete sequence, type strain: IAM 12609
AACTTGAGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGCCC
CGCAANNNNAGCGGCAGACGGGTGAGTAACGCGTGGGAACGTACCCTTTACTACGGAATAACGCAGGGAA
ACTTGTGCTAATACCGTATGTGCCCTTTGGGGGAAAGATTTATCGGTAAAGGATCGGCCCGCGTTGGATT
AGCTAGTTGGTGGGGTAAAGGCCTACCAAGGCGACGATCCATAGCTGGTCTGAGAGGATGATCAGCCACA
TTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGCAAGC
CTGATCCAGCCATGCCGCGTGAGTGATGAAGGCCCTAGGGTTGTAAAGCTCTTTCACCGGAGAAGATAAT
GACGGTATCCGGAGAAGAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCG
TTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGATCGATAAGTCAGGGGTGAAATCCCAGGGCTC
AACCCTGGAACTGCCTTTGATACTGTCGATCTGGAGTATGGAAGAGGTGAGTGGAATTCCGAGTGTAGAG
GTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGTCCATTACTGACGCTGAG
GTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATGTTAGC
CGTCGGGCAGTATACTGTTCGGTGGCGCACGTAACGCATTAAACATTCCGCCTGGGGAGTACGGTCGCAA
GATTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACG
CGCAGAACCTTACCAGCCCTTGACATGCCCGGCTACTTGCAGAGATGCAAGGTTCTTCGGGGACCGGGAC
ACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACC
CTCGCCCTTAGTTGCCAGCATTCAGTTGGGCACTCTAAGGGGACTGCCGGTGATAAGCCGAGAGGAAGGT
GGGGATGACGTCAAGTCCTCATGGCCCTTACGGGCTGGGCTACACACGTGCTACAATGGTGGTGACAGTG
GGCAGCGAGCACGCGAGTGTGAGCTAATCTCCAAAAGCCATCTCAGTTCGGATTGCACTCTGCAACTCGA
GTGCATGAAGTTGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTA
CACACCGCCCGTCACACCATGGGAGTTGGTTTTACCCGAAGGTAGTGCGCTAACCGCAAGGAGGCAGCTA
ACCACGGTAGGGTCAGCGACTGGGGTGAAGTCGTAACAAGGTAGCCGTAGGGGAACCTGCGGCTGGATCA
CCTCC

>NR_074106.1 Acidobacterium capsulatum ATCC 51196 16S ribosomal RNA, partial sequence
AGAGTTTGATCCTGGCTCAGAATCAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACAAGAAAGGGA
CTTCGGTCCTGAGTACAGTGGCGCACGGGTGAGTAACACGTGACTAACCTACCCTCGAGTGGGGAATAAC
TTCGGGAAACCGAGGCTAATACCGCATAATACCCACGGGTCAAAGGAGCAATTCGCTTGAGGAGGGGGTC
GCGGCCGATTAGCTAGTTGGCGGGGTAATGGCCCACCAAGGCAGTGATCGGTATCCGGCCTGAGAGGGCG
CACGGACACACTGGAACTGAAACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTGCGCAAT
GGGGGAAACCCTGACGCAGCAACGCCGCGTGGAGGATGAAGTCTCTTGGGACGTAAACTCCTTTCGATCG
GAACGATTATGACGGTACCGGAAGAAGAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAGG
GGGGCGAGCGTTGTTCGGAATTATTGGGCGTAAAGGGTGCGTAGGCGGTTCGGTAAGTTTGATGTGAAAT
CTTCGGGCTCAACTCGAAGTCTGCATCGAAAACTGCCGGGCTTGAGTGTGGGAGAGGTGAGTGGAATTTC
CGGTGTAGCGGTGAAATGCGTAGATATCGGAAGGAACACCTGTGGCGAAAGCGGCTCACTGGACCACAAC
TGACGCTGATGCACGAAAGCTAGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCTAGCCCTAAACGAT
GATCGCTTGGTGTGGCGGGTACCCAATCCCGTCGTGCCGTAGCTAACGCGTTAAGCGATCCGCCTGGGGA
GTACGGTCGCAAGGCTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAA
TTCGACGCAACGCGAAGAACCTTACCTGGGCTCGAAATGTAGTGGACCGGGGTAGAAATATCCCTTCCCC
GCAAGGGGCTGCTATATAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCC
GCAACGAGCGCAACCCTTATTGCCAGTTGCTACCATTTAGTTGAGCACTCTGGTGAGACCGCCTCGGATA
ACGGGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCTTTATGTCCAGGGCTACACACGTGCTACAA
TGGCCGGTACAAACCGCCGCAAACCCGCGAGGGGGAGCTAATCGGAAAAAGCCGGCCTCAGTTCGGATTG
TAGTCTGCAACTCGACTACATGAAGCTGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACG
TTCCCGGGCCTTGTACACACCGCCCGTCACATCACGAAAGTGGGTCGTACTAGAAGCGGGTGAGCCAACC
GTAAGGAGGCAGCCTTCCAAGGTGTGATTCATGATTGGGGTGAAGTCGTAACAAGGTAGCCGTAGGAGAA
CCTGCGGCTGGATCACCTCCTTT
"""


REFERENCE_META = {
    "NR_102783.2": {
        "label": "Bacillus_subtilis_168",
        "species": "Bacillus subtilis subsp. subtilis strain 168",
        "phylum": "Bacillota",
        "soil_context": "common soil and rhizosphere model bacterium",
        "note": "Gram-positive, spore-forming soil bacterium; useful classroom decomposer reference.",
    },
    "NR_115715.1": {
        "label": "Pseudomonas_fluorescens_CCM2115",
        "species": "Pseudomonas fluorescens strain CCM 2115",
        "phylum": "Pseudomonadota",
        "soil_context": "rhizosphere-associated soil bacterium",
        "note": "Common plant-root associated reference; useful contrast to Gram-positive taxa.",
    },
    "Y00411.1": {
        "label": "Streptomyces_coelicolor_rrnD",
        "species": "Streptomyces coelicolor",
        "phylum": "Actinomycetota",
        "soil_context": "filamentous soil actinomycete",
        "note": "Classic soil actinomycete; illustrates that soil microbes are not a single close group.",
    },
    "D14513.1": {
        "label": "Rhizobium_leguminosarum_IAM12609",
        "species": "Rhizobium leguminosarum type strain IAM 12609",
        "phylum": "Pseudomonadota",
        "soil_context": "root nodule and nitrogen-cycling context",
        "note": "Plant-associated nitrogen-cycle reference; sequence contains a few ambiguous N bases from the original record.",
    },
    "NR_074106.1": {
        "label": "Acidobacterium_capsulatum_ATCC51196",
        "species": "Acidobacterium capsulatum ATCC 51196",
        "phylum": "Acidobacteriota",
        "soil_context": "acidic soil and broad soil ecology reference",
        "note": "Soil-relevant reference from a major soil-associated phylum.",
    },
}


PHYLUM_COLORS = {
    "Bacillota": "#E69F00",
    "Pseudomonadota": "#56B4E9",
    "Actinomycetota": "#009E73",
    "Acidobacteriota": "#CC79A7",
    "Teaching query": "#222222",
}


def compact_sequence(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if line and not line.startswith(">")).upper()


def parse_fasta(text: str) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = None
    chunks: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header.split()[0]] = (header, "".join(chunks).upper())
            header = line[1:]
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        records[header.split()[0]] = (header, "".join(chunks).upper())
    return records


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def extract_marker(seq: str, max_len: int = 560) -> str:
    seq = seq.replace("U", "T").replace(" ", "").replace("\n", "").upper()
    anchor = "AGAGTTTGATCCTGGCTCAG"
    start = seq.find(anchor)
    if start < 0:
        start = max(seq.find("GAGTTTGATCCTGGCTCAG") - 1, 0)
    return seq[start : start + max_len]


def mutate(seq: str, positions: list[int]) -> str:
    bases = "ACGT"
    chars = list(seq)
    for pos in positions:
        if pos >= len(chars):
            continue
        current = chars[pos]
        if current not in bases:
            chars[pos] = "A"
        else:
            chars[pos] = bases[(bases.index(current) + 1) % len(bases)]
    return "".join(chars)


def global_align_pair(reference: str, sequence: str) -> tuple[str, str]:
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner(
        mode="global",
        match_score=2,
        mismatch_score=-1,
        open_gap_score=-5,
        extend_gap_score=-0.5,
    )
    alignment = aligner.align(reference, sequence)[0]
    ref_blocks, seq_blocks = alignment.aligned
    ref_parts: list[str] = []
    seq_parts: list[str] = []
    ref_pos = 0
    seq_pos = 0

    for (ref_start, ref_end), (seq_start, seq_end) in zip(ref_blocks, seq_blocks):
        if ref_start > ref_pos:
            ref_parts.append(reference[ref_pos:ref_start])
            seq_parts.append("-" * (ref_start - ref_pos))
        if seq_start > seq_pos:
            ref_parts.append("-" * (seq_start - seq_pos))
            seq_parts.append(sequence[seq_pos:seq_start])

        ref_parts.append(reference[ref_start:ref_end])
        seq_parts.append(sequence[seq_start:seq_end])
        ref_pos = int(ref_end)
        seq_pos = int(seq_end)

    if ref_pos < len(reference):
        ref_parts.append(reference[ref_pos:])
        seq_parts.append("-" * (len(reference) - ref_pos))
    if seq_pos < len(sequence):
        ref_parts.append("-" * (len(sequence) - seq_pos))
        seq_parts.append(sequence[seq_pos:])

    aligned_ref = "".join(ref_parts)
    aligned_seq = "".join(seq_parts)
    if len(aligned_ref) != len(aligned_seq):
        raise ValueError("PairwiseAligner returned unequal reconstructed alignment lengths")
    return aligned_ref, aligned_seq


def star_align_to_reference(windows: dict[str, str], reference_label: str) -> dict[str, str]:

    reference = windows[reference_label]
    ref_len = len(reference)
    base_by_label: dict[str, list[str]] = {}
    insertions_by_label: dict[str, dict[int, str]] = {}
    max_insertions: dict[int, int] = {i: 0 for i in range(ref_len + 1)}

    for label, sequence in windows.items():
        if label == reference_label:
            aligned_ref, aligned_seq = reference, sequence
        else:
            aligned_ref, aligned_seq = global_align_pair(reference, sequence)

        bases = ["-"] * ref_len
        insertions: dict[int, list[str]] = {}
        ref_pos = 0
        for ref_base, seq_base in zip(aligned_ref, aligned_seq):
            if ref_base == "-":
                insertions.setdefault(ref_pos, []).append(seq_base)
            else:
                if ref_pos < ref_len:
                    bases[ref_pos] = seq_base
                ref_pos += 1

        compact_insertions = {pos: "".join(chars) for pos, chars in insertions.items()}
        for pos, inserted in compact_insertions.items():
            max_insertions[pos] = max(max_insertions.get(pos, 0), len(inserted))
        base_by_label[label] = bases
        insertions_by_label[label] = compact_insertions

    aligned: dict[str, str] = {}
    for label in windows:
        chars: list[str] = []
        for pos in range(ref_len):
            inserted = insertions_by_label[label].get(pos, "")
            chars.append(inserted.ljust(max_insertions.get(pos, 0), "-"))
            chars.append(base_by_label[label][pos])
        terminal_insert = insertions_by_label[label].get(ref_len, "")
        chars.append(terminal_insert.ljust(max_insertions.get(ref_len, 0), "-"))
        aligned[label] = "".join(chars)

    lengths = {len(seq) for seq in aligned.values()}
    if len(lengths) != 1:
        raise ValueError(f"star alignment produced unequal lengths: {sorted(lengths)}")
    return aligned


def make_cache_files() -> dict[str, str]:
    raw_records = parse_fasta(RAW_REFERENCE_FASTA)

    ref_lines: list[str] = []
    metadata_rows: list[dict[str, str]] = []
    marker_by_accession: dict[str, str] = {}
    marker_by_label: dict[str, str] = {}

    for accession, (_header, seq) in raw_records.items():
        meta = REFERENCE_META[accession]
        label = meta["label"]
        ref_lines.append(
            f">{label} accession={accession} source=NCBI_Nucleotide role=reference species=\"{meta['species']}\""
        )
        ref_lines.append(wrap_fasta(seq))
        marker_by_accession[accession] = extract_marker(seq)
        marker_by_label[label] = marker_by_accession[accession]
        metadata_rows.append(
            {
                "label": label,
                "role": "reference",
                "species_or_query": meta["species"],
                "accession": accession,
                "source_database": "NCBI Nucleotide",
                "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                "retrieved_date": RETRIEVED_DATE,
                "phylum": meta["phylum"],
                "color": PHYLUM_COLORS[meta["phylum"]],
                "soil_context": meta["soil_context"],
                "note": meta["note"],
            }
        )

    query_records = [
        {
            "label": "Soil_ASV_A",
            "source_accession": "NR_102783.2",
            "species_or_query": "unknown soil ASV A",
            "phylum": "Teaching query",
            "soil_context": "class teaching query from a rhizosphere-style sample",
            "note": "Synthetic classroom read derived from the Bacillus subtilis 16S marker window with a few substitutions.",
            "sequence": mutate(marker_by_accession["NR_102783.2"][:520], [43, 141, 309, 471]),
        },
        {
            "label": "Soil_ASV_B",
            "source_accession": "D14513.1",
            "species_or_query": "unknown soil ASV B",
            "phylum": "Teaching query",
            "soil_context": "class teaching query from a root-associated soil sample",
            "note": "Synthetic classroom read derived from the Rhizobium leguminosarum 16S marker window with a few substitutions.",
            "sequence": mutate(marker_by_accession["D14513.1"][:520], [58, 220, 391, 505]),
        },
    ]

    def distance_stats(seq_a: str, seq_b: str) -> tuple[int, int, float]:
        n = min(len(seq_a), len(seq_b))
        compared = 0
        differences = 0
        for left, right in zip(seq_a[:n], seq_b[:n]):
            if left in "N-" or right in "N-":
                continue
            compared += 1
            if left != right:
                differences += 1
        fraction = differences / compared if compared else math.nan
        return compared, differences, fraction

    cached_hit_rows: list[dict[str, str]] = []
    window_by_label = dict(marker_by_label)
    window_by_label.update({query["label"]: query["sequence"] for query in query_records})
    aligned_for_hits = star_align_to_reference(window_by_label, "Bacillus_subtilis_168")
    for query in query_records:
        ranked_hits: list[dict[str, str]] = []
        for accession, meta in REFERENCE_META.items():
            compared, differences, fraction = distance_stats(
                aligned_for_hits[query["label"]],
                aligned_for_hits[meta["label"]],
            )
            ranked_hits.append(
                {
                    "query_label": query["label"],
                    "reference_label": meta["label"],
                    "reference_accession": accession,
                    "reference_species": meta["species"],
                    "compared_bases": str(compared),
                    "differences": str(differences),
                    "fraction_different": f"{fraction:.6f}",
                    "percent_identity_teaching_window": f"{100 * (1 - fraction):.3f}",
                    "source": "precomputed class-safe teaching hit table",
                }
            )
        ranked_hits.sort(key=lambda row: float(row["fraction_different"]))
        for rank, row in enumerate(ranked_hits, start=1):
            row["rank"] = str(rank)
            row["interpretation"] = (
                "closest cached reference"
                if rank == 1
                else "lower-ranked cached reference"
            )
            cached_hit_rows.append(row)

    query_lines: list[str] = []
    for query in query_records:
        query_lines.append(
            f">{query['label']} source_accession={query['source_accession']} source=teaching_cache role=query"
        )
        query_lines.append(wrap_fasta(query["sequence"]))
        metadata_rows.append(
            {
                "label": query["label"],
                "role": "query",
                "species_or_query": query["species_or_query"],
                "accession": query["source_accession"],
                "source_database": "Teaching cache derived from NCBI reference",
                "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{query['source_accession']}",
                "retrieved_date": RETRIEVED_DATE,
                "phylum": query["phylum"],
                "color": PHYLUM_COLORS["Teaching query"],
                "soil_context": query["soil_context"],
                "note": query["note"],
            }
        )

    metadata_fields = [
        "label",
        "role",
        "species_or_query",
        "accession",
        "source_database",
        "source_url",
        "retrieved_date",
        "phylum",
        "color",
        "soil_context",
        "note",
    ]
    metadata_csv = rows_to_csv(metadata_rows, metadata_fields)

    cached_hits_fields = [
        "query_label",
        "rank",
        "reference_label",
        "reference_accession",
        "reference_species",
        "compared_bases",
        "differences",
        "fraction_different",
        "percent_identity_teaching_window",
        "source",
        "interpretation",
    ]
    cached_hits_csv = rows_to_csv(cached_hit_rows, cached_hits_fields)
    cached_blast_xml = cached_hits_to_blast_xml(cached_hit_rows)

    abundance_rows = [
        {"sample_id": "Rhizosphere_A", "Soil_ASV_A": "128", "Soil_ASV_B": "34", "note": "Bacillus-like ASV is more abundant in this toy sample."},
        {"sample_id": "Compost_B", "Soil_ASV_A": "76", "Soil_ASV_B": "93", "note": "Both ASVs are detectable in the compost-style toy sample."},
        {"sample_id": "Root_Nodule_C", "Soil_ASV_A": "21", "Soil_ASV_B": "156", "note": "Rhizobium-like ASV is more abundant in this toy sample."},
    ]
    abundance_csv = rows_to_csv(abundance_rows, ["sample_id", "Soil_ASV_A", "Soil_ASV_B", "note"])

    manifest = {
        "title": "Class-safe soil 16S phylogeny pilot cache",
        "retrieved_date": RETRIEVED_DATE,
        "reference_count": 5,
        "query_count": 2,
        "cached_hit_table": "pilot_16s_cached_hits.csv",
        "cached_blast_xml": "pilot_16s_cached_blast.xml",
        "default_mode": "Use these cached files; do not require live BLAST or Entrez during class.",
        "references": [
            {
                "label": meta["label"],
                "accession": accession,
                "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
            }
            for accession, meta in REFERENCE_META.items()
        ],
        "teaching_queries": [
            {
                "label": query["label"],
                "derived_from": query["source_accession"],
                "caution": "Synthetic classroom query; not a new environmental isolate.",
            }
            for query in query_records
        ],
    }

    readme = textwrap.dedent(
        f"""
        # Soil 16S Class Cache

        This folder is the class-safe data cache for `soil_microbiome_16s_class_safe_colab.ipynb`.

        The default classroom workflow should load these files instead of calling live BLAST, Entrez, SILVA, or other web services during class.

        This cache supports a 16S marker/metabarcoding teaching workflow, not a shotgun metagenomics workflow.

        One-command Colab pattern after this folder is pushed to GitHub:

        ```python
        CACHE_BASE_URL = "https://raw.githubusercontent.com/<org>/<repo>/main/soil_16s_class_cache"
        ```

        Files:

        - `pilot_16s_references.fasta`: five NCBI 16S reference sequences.
        - `pilot_16s_query_reads.fasta`: two synthetic teaching ASV/read sequences derived from the cached references.
        - `pilot_16s_metadata.csv`: source, accession, phylum, soil context, and interpretation notes.
        - `pilot_16s_cached_hits.csv`: precomputed closest-reference table for the teaching queries.
        - `pilot_16s_cached_blast.xml`: cached BLAST-like XML for teaching ranked-hit parsing without live BLAST.
        - `pilot_16s_abundance_table.csv`: toy microbiome-style counts.
        - `pilot_16s_manifest.json`: provenance and class-use notes.
        - `cache_validation_report.json`: last local validation result from the builder.
        - `notebook_execution_report.json`: result from executing the generated notebook with the verifier.
        - `VISUAL_QA_TUFTE.md`: visualization checklist used for the student-facing figures.
        - `COLAB_ONE_CELL_LOADER.py`: copy-paste loader for any other Colab notebook.
        - `PUBLISH_TO_GITHUB.md`: publish and pre-class verification checklist.

        Retrieved date recorded in the manifest: {RETRIEVED_DATE}.
        """
    ).strip() + "\n"

    loader_py = '''"""One-cell loader for the class-safe soil 16S cache.

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
'''

    publish_md = """# Publish The Soil 16S Class Cache

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
"""

    visual_qa_md = """# Visual QA: Tufte Checklist

This checklist applies to `soil_microbiome_16s_class_safe_colab.ipynb`.

## Palette

- Categories use Okabe-Ito colors.
- Continuous distance values use `viridis`.
- No decorative gradients, glow, or heavy card styling.

## Figure Checks

- Workflow map: labels sit above the marks; connectors are thin gray lines.
- Abundance plot: stacked bars encode the toy count table; no 3D effects.
- Alignment view: base colors encode A/C/G/T/N/gap and include a compact legend; variable columns are small ticks.
- Distance matrix: colorbar says exactly what the values mean.
- Tree plots: branch length axis remains visible; unnecessary plot borders are removed.

## Eraser Test

Every visible element should either encode data, orient the student, or label the data. Remove decorations that do none of these jobs.

## Collision Test

Before teaching, run the notebook in Colab at normal width and inspect:

- axis labels,
- tip labels,
- workflow labels,
- table wrapping,
- heatmap tick labels.

No label should overlap another label or hide a data mark.

## Scientific Integrity

- The distance heatmap uses a zero baseline and a monotonic sequential colormap.
- The notebook states that the tree is built from one 16S marker window.
- The report sentence says closest reference, not exact species proof.
"""

    return {
        "pilot_16s_references.fasta": "\n".join(ref_lines).strip() + "\n",
        "pilot_16s_query_reads.fasta": "\n".join(query_lines).strip() + "\n",
        "pilot_16s_metadata.csv": metadata_csv,
        "pilot_16s_cached_hits.csv": cached_hits_csv,
        "pilot_16s_cached_blast.xml": cached_blast_xml,
        "pilot_16s_abundance_table.csv": abundance_csv,
        "pilot_16s_manifest.json": json.dumps(manifest, indent=2) + "\n",
        "README.md": readme,
        "COLAB_ONE_CELL_LOADER.py": loader_py,
        "PUBLISH_TO_GITHUB.md": publish_md,
        "VISUAL_QA_TUFTE.md": visual_qa_md,
    }


def rows_to_csv(rows: list[dict[str, str]], fields: list[str]) -> str:
    from io import StringIO

    buffer = StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def cached_hits_to_blast_xml(rows: list[dict[str, str]]) -> str:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["query_label"], []).append(row)

    lines = [
        '<?xml version="1.0"?>',
        "<BlastOutput>",
        "  <BlastOutput_program>blastn</BlastOutput_program>",
        "  <BlastOutput_version>cached-teaching-blast-1.0</BlastOutput_version>",
        "  <BlastOutput_reference>Class-safe teaching XML generated from precomputed aligned-window hits; not a live NCBI BLAST run.</BlastOutput_reference>",
        "  <BlastOutput_db>soil_16s_class_cache</BlastOutput_db>",
        "  <BlastOutput_query-ID>soil_16s_teaching_queries</BlastOutput_query-ID>",
        "  <BlastOutput_query-def>soil_16s_teaching_queries</BlastOutput_query-def>",
        "  <BlastOutput_param><Parameters><Parameters_matrix>identity</Parameters_matrix></Parameters></BlastOutput_param>",
        "  <BlastOutput_iterations>",
    ]
    for iter_num, (query_label, query_rows) in enumerate(grouped.items(), start=1):
        query_rows = sorted(query_rows, key=lambda item: int(item["rank"]))
        max_compared = max(int(row["compared_bases"]) for row in query_rows)
        lines.extend(
            [
                "    <Iteration>",
                f"      <Iteration_iter-num>{iter_num}</Iteration_iter-num>",
                f"      <Iteration_query-ID>{xml_escape(query_label)}</Iteration_query-ID>",
                f"      <Iteration_query-def>{xml_escape(query_label)}</Iteration_query-def>",
                f"      <Iteration_query-len>{max_compared}</Iteration_query-len>",
                "      <Iteration_hits>",
            ]
        )
        for row in query_rows:
            compared = int(row["compared_bases"])
            differences = int(row["differences"])
            identity_count = compared - differences
            bit_score = round(2 * identity_count - differences, 3)
            evalue = "1e-120" if row["rank"] == "1" else f"1e-{max(10, 90 - 5 * int(row['rank']))}"
            lines.extend(
                [
                    "        <Hit>",
                    f"          <Hit_num>{row['rank']}</Hit_num>",
                    f"          <Hit_id>{xml_escape(row['reference_label'])}</Hit_id>",
                    f"          <Hit_def>{xml_escape(row['reference_species'])}</Hit_def>",
                    f"          <Hit_accession>{xml_escape(row['reference_accession'])}</Hit_accession>",
                    f"          <Hit_len>{compared}</Hit_len>",
                    "          <Hit_hsps>",
                    "            <Hsp>",
                    f"              <Hsp_num>1</Hsp_num>",
                    f"              <Hsp_bit-score>{bit_score}</Hsp_bit-score>",
                    f"              <Hsp_evalue>{evalue}</Hsp_evalue>",
                    f"              <Hsp_query-from>1</Hsp_query-from>",
                    f"              <Hsp_query-to>{compared}</Hsp_query-to>",
                    f"              <Hsp_hit-from>1</Hsp_hit-from>",
                    f"              <Hsp_hit-to>{compared}</Hsp_hit-to>",
                    f"              <Hsp_identity>{identity_count}</Hsp_identity>",
                    f"              <Hsp_gaps>0</Hsp_gaps>",
                    f"              <Hsp_align-len>{compared}</Hsp_align-len>",
                    "            </Hsp>",
                    "          </Hit_hsps>",
                    "        </Hit>",
                ]
            )
        lines.extend(["      </Iteration_hits>", "    </Iteration>"])
    lines.extend(["  </BlastOutput_iterations>", "</BlastOutput>", ""])
    return "\n".join(lines)


def md(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": textwrap.dedent(source).strip() + "\n"}


def code(source: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": textwrap.dedent(source).strip() + "\n",
    }


def make_notebook(cache_files: dict[str, str]) -> dict:
    cache_literal = json.dumps(
        {name: content for name, content in cache_files.items() if name.startswith("pilot_")},
        indent=2,
    ).replace("\n", "\n            ")

    cells = [
        md(
            """
            # Soil 16S Phylogeny Pilot

            ## Species finding and relatedness from a class-safe cache

            Run the cells from top to bottom. When you see a **Think** prompt, pause and write one sentence.

            Today you will:

            1. load cached 16S rRNA reference sequences from soil-relevant bacteria,
            2. compare two unknown soil ASV/read sequences with those references,
            3. inspect a small aligned marker window,
            4. turn sequence differences into a distance matrix,
            5. build UPGMA and neighbor-joining trees,
            6. report a careful closest-reference claim.

            A tree is a hypothesis from evidence. In this notebook the evidence is one short 16S marker window, so the final claim must stay cautious.

            This is a browser-only 16S marker/metabarcoding pilot. It is not a shotgun metagenomics pipeline and it does not prove exact species identity.
            """
        ),
        code(
            """
            #@title Class controls { display-mode: "form" }
            USE_GITHUB_CACHE = False #@param {type:"boolean"}
            CACHE_BASE_URL = "" #@param {type:"string"}
            QUERY_TO_REPORT = "Soil_ASV_A" #@param ["Soil_ASV_A", "Soil_ASV_B"]
            TREE_METHOD_TO_SHOW = "Compare UPGMA and neighbor joining" #@param ["UPGMA", "Neighbor joining", "Compare UPGMA and neighbor joining"]
            MARKER_WINDOW_BASES = 520 #@param {type:"slider", min:200, max:560, step:20}
            ALIGNMENT_START = 130 #@param {type:"slider", min:0, max:420, step:10}
            ALIGNMENT_WIDTH = 70 #@param {type:"slider", min:30, max:100, step:10}
            print("Controls set. The default path uses the embedded cache, so class runs do not depend on live BLAST.")
            """
        ),
        code(
            f"""
            #@title Install and import notebook dependencies {{ display-mode: "form" }}
            import importlib.util
            import subprocess
            import sys

            def ensure(package, import_name=None):
                import_name = import_name or package
                if importlib.util.find_spec(import_name) is None:
                    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", package])

            for package, import_name in [
                ("biopython", "Bio"),
                ("pandas", "pandas"),
                ("matplotlib", "matplotlib"),
                ("seaborn", "seaborn"),
            ]:
                ensure(package, import_name)

            import json
            import math
            import urllib.request
            import xml.etree.ElementTree as ET
            from io import StringIO
            from pathlib import Path

            import numpy as np
            import pandas as pd
            import matplotlib.pyplot as plt
            import seaborn as sns
            from IPython.display import HTML, display

            from Bio import Phylo, SeqIO
            from Bio.Align import PairwiseAligner
            from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

            EMBEDDED_CACHE = {cache_literal}

            OKABE_ITO = {{
                "orange": "#E69F00",
                "sky_blue": "#56B4E9",
                "bluish_green": "#009E73",
                "yellow": "#F0E442",
                "blue": "#0072B2",
                "vermillion": "#D55E00",
                "reddish_purple": "#CC79A7",
                "black": "#222222",
                "gray": "#DDDDDD",
            }}

            BASE_COLORS = {{
                "A": OKABE_ITO["sky_blue"],
                "C": OKABE_ITO["bluish_green"],
                "G": OKABE_ITO["orange"],
                "T": OKABE_ITO["reddish_purple"],
                "N": "#BDBDBD",
                "-": "#E6E6E6",
            }}

            plt.rcParams.update({{
                "figure.facecolor": "white",
                "axes.facecolor": "white",
                "axes.edgecolor": "#444444",
                "axes.labelcolor": "#222222",
                "text.color": "#222222",
                "xtick.color": "#222222",
                "ytick.color": "#222222",
                "font.family": "DejaVu Sans",
                "font.size": 10,
                "axes.titleweight": "regular",
                "savefig.facecolor": "white",
            }})
            sns.set_theme(style="white")

            print("Ready. Dependencies are available and the embedded class cache is loaded in memory.")
            """
        ),
        md(
            """
            ## 1. Why cached data?

            What can fail during a live class? A public database can be slow, a web API can throttle, or a student can lose time debugging a connection.

            So this pilot uses a pre-cached reference set. The biology is still real enough for the lesson: the five reference records are NCBI 16S rRNA sequences from soil-relevant bacteria. The two query reads are clearly marked teaching reads derived from that cache.

            **Think:** Why is cached data better for the first teaching run than live BLAST?
            """
        ),
        code(
            """
            #@title Load cache from GitHub or embedded fallback
            def load_cache_file(filename):
                if USE_GITHUB_CACHE and CACHE_BASE_URL.strip():
                    url = CACHE_BASE_URL.rstrip("/") + "/" + filename
                    try:
                        with urllib.request.urlopen(url, timeout=20) as handle:
                            text = handle.read().decode("utf-8")
                        print(f"Loaded {filename} from GitHub cache.")
                        return text
                    except Exception as exc:
                        print(f"GitHub cache failed for {filename}: {exc}")
                        print("Using embedded fallback copy.")
                return EMBEDDED_CACHE[filename]

            ref_fasta = load_cache_file("pilot_16s_references.fasta")
            query_fasta = load_cache_file("pilot_16s_query_reads.fasta")
            metadata_csv = load_cache_file("pilot_16s_metadata.csv")
            cached_hits_csv = load_cache_file("pilot_16s_cached_hits.csv")
            cached_blast_xml = load_cache_file("pilot_16s_cached_blast.xml")
            abundance_csv = load_cache_file("pilot_16s_abundance_table.csv")
            manifest_json = load_cache_file("pilot_16s_manifest.json")

            references = list(SeqIO.parse(StringIO(ref_fasta), "fasta"))
            queries = list(SeqIO.parse(StringIO(query_fasta), "fasta"))
            all_records = references + queries
            metadata = pd.read_csv(StringIO(metadata_csv))
            cached_hits = pd.read_csv(StringIO(cached_hits_csv))
            abundance = pd.read_csv(StringIO(abundance_csv))
            manifest = json.loads(manifest_json)

            def parse_cached_blast_xml(xml_text):
                root = ET.fromstring(xml_text)
                rows = []
                for iteration in root.findall(".//Iteration"):
                    query = iteration.findtext("Iteration_query-def")
                    for hit in iteration.findall("./Iteration_hits/Hit"):
                        hsp = hit.find("./Hit_hsps/Hsp")
                        align_len = int(hsp.findtext("Hsp_align-len"))
                        identity = int(hsp.findtext("Hsp_identity"))
                        rows.append({
                            "query_label": query,
                            "rank": int(hit.findtext("Hit_num")),
                            "reference_label": hit.findtext("Hit_id"),
                            "reference_accession": hit.findtext("Hit_accession"),
                            "reference_species": hit.findtext("Hit_def"),
                            "compared_bases": align_len,
                            "differences": align_len - identity,
                            "percent_identity_teaching_window": 100 * identity / align_len,
                            "bit_score": float(hsp.findtext("Hsp_bit-score")),
                            "teaching_e_value": hsp.findtext("Hsp_evalue"),
                        })
                return pd.DataFrame(rows)

            cached_blast_hits = parse_cached_blast_xml(cached_blast_xml)

            print(f"Reference sequences: {len(references)}")
            print(f"Query reads: {len(queries)}")
            print(f"Cached BLAST-like XML hits: {len(cached_blast_hits)}")
            print("Cache retrieved date:", manifest["retrieved_date"])
            """
        ),
        code(
            """
            #@title Show the teaching reference set
            def wrapped_table(df, columns):
                rows = []
                for _, row in df[columns].iterrows():
                    cells = "".join(
                        f"<td style='padding:8px 10px; border-top:1px solid #e5e5e5; vertical-align:top;'>{row[col]}</td>"
                        for col in columns
                    )
                    rows.append(f"<tr>{cells}</tr>")
                header = "".join(
                    f"<th style='text-align:left; padding:7px 10px; border-bottom:1px solid #777;'>{col}</th>"
                    for col in columns
                )
                html = f'''
                <div style='font-family: system-ui, Segoe UI, sans-serif; max-width: 980px;'>
                  <table style='border-collapse: collapse; table-layout: fixed; width: 100%; font-size: 13px; line-height: 1.35;'>
                    <thead><tr>{header}</tr></thead>
                    <tbody>{''.join(rows)}</tbody>
                  </table>
                </div>
                '''
                display(HTML(html))

            display_cols = ["label", "role", "species_or_query", "phylum", "soil_context", "note"]
            wrapped_table(metadata, display_cols)
            """
        ),
        md(
            """
            ## 2. The class-safe workflow

            What should happen when a student presses **Run all**?

            The notebook should load known files, make the same comparison for everyone, and produce interpretable figures. Live BLAST and live database searches are useful later, but they are not allowed to be the first-class path because they can fail during class.

            The workflow is:

            **cached references + cached query reads -> marker window -> alignment view -> distance matrix -> UPGMA/NJ trees -> cautious report**
            """
        ),
        code(
            """
            #@title Workflow map
            steps = [
                ("cache", "cached FASTA + metadata"),
                ("window", "shared 16S marker window"),
                ("align", "aligned position view"),
                ("distance", "pairwise distances"),
                ("tree", "UPGMA + NJ trees"),
                ("claim", "closest-reference report"),
            ]

            fig, ax = plt.subplots(figsize=(11, 2.2))
            ax.set_xlim(-0.5, len(steps) - 0.5)
            ax.set_ylim(-0.5, 1.05)
            ax.axis("off")

            for i, (short, label) in enumerate(steps):
                ax.scatter(i, 0.35, s=160, color=OKABE_ITO["sky_blue"], edgecolor="#222222", linewidth=0.8, zorder=3)
                ax.text(i, 0.72, label, ha="center", va="bottom", fontsize=9, wrap=True)
                ax.text(i, 0.35, str(i + 1), ha="center", va="center", fontsize=9, color="white", weight="bold")
                if i < len(steps) - 1:
                    ax.plot([i + 0.15, i + 0.85], [0.35, 0.35], color="#999999", linewidth=1.0, zorder=1)

            ax.set_title("Class-safe path: no live web service is required during the lesson", loc="left", fontsize=12)
            plt.tight_layout()
            plt.show()
            """
        ),
        md(
            """
            ## 3. Where the references come from

            In a real soil microbiome project, a 16S read can be compared against databases such as NCBI BLAST/GenBank, SILVA, RDP, and GTDB.

            For this pilot, we already did the retrieval step before class and cached a small reference set. That makes the first run reliable. Later, the class can replace the teaching cache with references selected from live searches.
            """
        ),
        code(
            """
            #@title Reference database map
            database_rows = pd.DataFrame([
                {
                    "resource": "NCBI BLAST / GenBank",
                    "helps with": "finding close public sequence records and accessions",
                    "class-safe use": "cache selected FASTA and source links before class",
                },
                {
                    "resource": "SILVA",
                    "helps with": "curated rRNA reference alignment and taxonomy context",
                    "class-safe use": "download or export selected references before class",
                },
                {
                    "resource": "RDP",
                    "helps with": "ribosomal RNA taxonomy and classifier-style teaching comparisons",
                    "class-safe use": "cache selected reference records or classifier output",
                },
                {
                    "resource": "GTDB",
                    "helps with": "modern genome-based bacterial taxonomy context",
                    "class-safe use": "use as taxonomy background, not as a live dependency",
                },
            ])
            wrapped_table(database_rows, ["resource", "helps with", "class-safe use"])
            """
        ),
        md(
            """
            ## 4. Cached closest-hit table

            What would students normally wait for from BLAST?

            They would wait for a ranked hit table. For class, we cache that idea too. This table is computed from the same teaching marker window and lets students see the species-finding result before they build a tree.

            The E-values and bit scores shown below come from the cached BLAST-like teaching XML. They are included so students can learn how ranked-hit evidence is read, but they are not live NCBI BLAST statistics.
            """
        ),
        code(
            """
            #@title Cached closest-hit table
            xml_top = (
                cached_blast_hits[cached_blast_hits["rank"] == 1]
                .set_index("query_label")["reference_label"]
                .to_dict()
            )
            csv_top = (
                cached_hits[cached_hits["rank"] == 1]
                .set_index("query_label")["reference_label"]
                .to_dict()
            )
            assert xml_top == csv_top, f"Cached XML and CSV disagree: {xml_top} vs {csv_top}"

            hit_view = (
                cached_blast_hits[cached_blast_hits["rank"] <= 3]
                .merge(
                    cached_hits[["query_label", "rank", "fraction_different", "interpretation"]],
                    on=["query_label", "rank"],
                    how="left",
                )
                .copy()
            )
            hit_view["percent_identity_teaching_window"] = hit_view["percent_identity_teaching_window"].map(lambda x: f"{float(x):.2f}")
            hit_view["fraction_different"] = hit_view["fraction_different"].map(lambda x: f"{float(x):.4f}")
            hit_view["bit_score"] = hit_view["bit_score"].map(lambda x: f"{float(x):.1f}")
            wrapped_table(
                hit_view[
                    [
                        "query_label",
                        "rank",
                        "reference_label",
                        "reference_accession",
                        "percent_identity_teaching_window",
                        "bit_score",
                        "teaching_e_value",
                        "interpretation",
                    ]
                ],
                [
                    "query_label",
                    "rank",
                    "reference_label",
                    "reference_accession",
                    "percent_identity_teaching_window",
                    "bit_score",
                    "teaching_e_value",
                    "interpretation",
                ],
            )
            print("Cached BLAST-like XML agrees with the cached hit CSV for the top hit of each query.")
            print("E-values and bit scores here are cached teaching values, not fresh live-BLAST statistics.")
            """
        ),
        md(
            """
            ## 5. From soil sample to species-finding question

            A microbiome table can tell us that a read or ASV exists in a sample. It does not automatically tell us what species it is.

            The next question is narrower: which cached reference sequence is each query most similar to?

            We answer that with sequence comparison first, then we use a tree to report relatedness.
            """
        ),
        code(
            """
            #@title Tiny microbiome-style count table
            fig, ax = plt.subplots(figsize=(8.5, 3.2))
            counts = abundance.set_index("sample_id")[["Soil_ASV_A", "Soil_ASV_B"]]
            counts.plot(
                kind="bar",
                stacked=True,
                color=[OKABE_ITO["blue"], OKABE_ITO["orange"]],
                edgecolor="white",
                linewidth=0.7,
                ax=ax,
            )
            ax.set_ylabel("teaching read count")
            ax.set_xlabel("")
            ax.set_title("Toy soil samples: two ASVs to identify", loc="left", fontsize=12)
            ax.legend(frameon=False, title="")
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="y", color="#eeeeee", linewidth=0.8)
            plt.xticks(rotation=0, ha="center")
            plt.tight_layout()
            plt.show()
            wrapped_table(abundance, ["sample_id", "Soil_ASV_A", "Soil_ASV_B", "note"])
            """
        ),
        md(
            """
            ## 6. Make a comparable 16S marker window

            How can we compare sequences fairly? First, we choose the same marker region.

            This pilot trims each reference to a shared 16S starting anchor and then uses Biopython's `PairwiseAligner` to make a small star alignment around one reference coordinate. That keeps the class path simple, but still makes the distance step depend on aligned columns rather than raw string positions.

            In the later project notebook, this shortcut can be replaced by MAFFT on the team's cleaned reads and selected references.
            """
        ),
        code(
            """
            #@title Extract the class marker window
            def clean_sequence(record):
                return str(record.seq).upper().replace("U", "T").replace(" ", "")

            def marker_window(seq, max_bases):
                anchor = "AGAGTTTGATCCTGGCTCAG"
                start = seq.find(anchor)
                if start < 0:
                    backup = seq.find("GAGTTTGATCCTGGCTCAG")
                    start = max(backup - 1, 0) if backup >= 0 else 0
                return seq[start : start + max_bases]

            windows = {}
            for record in all_records:
                seq = clean_sequence(record)
                if record.id.startswith("Soil_ASV"):
                    windows[record.id] = seq[:MARKER_WINDOW_BASES]
                else:
                    windows[record.id] = marker_window(seq, MARKER_WINDOW_BASES)

            def global_align_pair(reference, sequence):
                aligner = PairwiseAligner(
                    mode="global",
                    match_score=2,
                    mismatch_score=-1,
                    open_gap_score=-5,
                    extend_gap_score=-0.5,
                )
                alignment = aligner.align(reference, sequence)[0]
                ref_blocks, seq_blocks = alignment.aligned
                ref_parts = []
                seq_parts = []
                ref_pos = 0
                seq_pos = 0

                for (ref_start, ref_end), (seq_start, seq_end) in zip(ref_blocks, seq_blocks):
                    if ref_start > ref_pos:
                        ref_parts.append(reference[ref_pos:ref_start])
                        seq_parts.append("-" * (ref_start - ref_pos))
                    if seq_start > seq_pos:
                        ref_parts.append("-" * (seq_start - seq_pos))
                        seq_parts.append(sequence[seq_pos:seq_start])

                    ref_parts.append(reference[ref_start:ref_end])
                    seq_parts.append(sequence[seq_start:seq_end])
                    ref_pos = int(ref_end)
                    seq_pos = int(seq_end)

                if ref_pos < len(reference):
                    ref_parts.append(reference[ref_pos:])
                    seq_parts.append("-" * (len(reference) - ref_pos))
                if seq_pos < len(sequence):
                    ref_parts.append("-" * (len(sequence) - seq_pos))
                    seq_parts.append(sequence[seq_pos:])

                aligned_ref = "".join(ref_parts)
                aligned_seq = "".join(seq_parts)
                if len(aligned_ref) != len(aligned_seq):
                    raise ValueError("PairwiseAligner returned unequal reconstructed alignment lengths")
                return aligned_ref, aligned_seq


            def star_align_to_reference(windows, reference_label):
                reference = windows[reference_label]
                ref_len = len(reference)
                base_by_label = {}
                insertions_by_label = {}
                max_insertions = {i: 0 for i in range(ref_len + 1)}

                for label, sequence in windows.items():
                    if label == reference_label:
                        aligned_ref, aligned_seq = reference, sequence
                    else:
                        aligned_ref, aligned_seq = global_align_pair(reference, sequence)

                    bases = ["-"] * ref_len
                    insertions = {}
                    ref_pos = 0
                    for ref_base, seq_base in zip(aligned_ref, aligned_seq):
                        if ref_base == "-":
                            insertions.setdefault(ref_pos, []).append(seq_base)
                        else:
                            if ref_pos < ref_len:
                                bases[ref_pos] = seq_base
                            ref_pos += 1

                    compact_insertions = {pos: "".join(chars) for pos, chars in insertions.items()}
                    for pos, inserted in compact_insertions.items():
                        max_insertions[pos] = max(max_insertions.get(pos, 0), len(inserted))
                    base_by_label[label] = bases
                    insertions_by_label[label] = compact_insertions

                aligned = {}
                for label in windows:
                    chars = []
                    for pos in range(ref_len):
                        inserted = insertions_by_label[label].get(pos, "")
                        chars.append(inserted.ljust(max_insertions.get(pos, 0), "-"))
                        chars.append(base_by_label[label][pos])
                    terminal_insert = insertions_by_label[label].get(ref_len, "")
                    chars.append(terminal_insert.ljust(max_insertions.get(ref_len, 0), "-"))
                    aligned[label] = "".join(chars)
                return aligned

            aligned_windows = star_align_to_reference(windows, "Bacillus_subtilis_168")

            window_lengths = pd.DataFrame({
                "raw_marker_bases": {label: len(seq) for label, seq in windows.items()},
                "aligned_columns": {label: len(seq) for label, seq in aligned_windows.items()},
            })
            display(window_lengths)
            if window_lengths["raw_marker_bases"].min() < MARKER_WINDOW_BASES:
                print("Some records are shorter than the selected window; distances use the available bases.")
            """
        ),
        code(
            """
            #@title Visualize a small alignment window
            from matplotlib.colors import ListedColormap
            from matplotlib.patches import Patch

            def plot_alignment_window(windows, start=0, width=70):
                labels = list(windows)
                end = min(start + width, min(len(seq) for seq in windows.values()))
                bases = ["A", "C", "G", "T", "N", "-"]
                base_to_int = {base: i for i, base in enumerate(bases)}
                matrix = np.array([
                    [base_to_int.get(base, base_to_int["N"]) for base in windows[label][start:end]]
                    for label in labels
                ])

                fig, ax = plt.subplots(figsize=(11, 4.1))
                cmap = ListedColormap([BASE_COLORS[base] for base in bases])
                ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=len(bases)-1)
                ax.set_yticks(range(len(labels)))
                ax.set_yticklabels(labels)
                ax.set_xticks(range(0, end - start, 10))
                ax.set_xticklabels([str(start + x) for x in range(0, end - start, 10)])
                ax.set_xlabel("aligned marker-window column")
                ax.set_title("Aligned teaching window: conserved columns stay quiet; variable columns carry the signal", loc="left", fontsize=12)
                ax.tick_params(length=0)
                for spine in ax.spines.values():
                    spine.set_visible(False)
                legend_handles = [Patch(facecolor=BASE_COLORS[base], edgecolor="none", label=base) for base in bases]
                ax.legend(
                    handles=legend_handles,
                    loc="upper center",
                    bbox_to_anchor=(0.5, -0.2),
                    ncol=len(bases),
                    frameon=False,
                    handlelength=1.0,
                    columnspacing=1.0,
                )

                for x in range(matrix.shape[1]):
                    column = [windows[label][start + x] for label in labels if start + x < len(windows[label])]
                    observed = {base for base in column if base in "ACGT"}
                    if len(observed) > 1:
                        ax.plot([x, x], [-0.48, -0.25], color="#222222", linewidth=0.7, clip_on=False)

                plt.tight_layout()
                plt.show()

            plot_alignment_window(aligned_windows, ALIGNMENT_START, ALIGNMENT_WIDTH)
            """
        ),
        md(
            """
            ## 7. Compute sequence distances

            Once the sequences are comparable, what is the simplest number we can calculate?

            For each pair, we count how many aligned positions differ. Dividing by the number of compared positions gives a fraction different. Small values mean the two marker sequences are similar. Larger values mean they are more different.
            """
        ),
        code(
            """
            #@title Pairwise distances and closest references
            def fraction_different(seq_a, seq_b):
                n = min(len(seq_a), len(seq_b))
                a = seq_a[:n]
                b = seq_b[:n]
                compared = 0
                differences = 0
                for left, right in zip(a, b):
                    if left in "N-" or right in "N-":
                        continue
                    compared += 1
                    if left != right:
                        differences += 1
                return differences / compared if compared else math.nan

            labels = list(aligned_windows)
            dist = pd.DataFrame(index=labels, columns=labels, dtype=float)
            for left in labels:
                for right in labels:
                    dist.loc[left, right] = fraction_different(aligned_windows[left], aligned_windows[right])

            query_labels = [record.id for record in queries]
            reference_labels = [record.id for record in references]
            tree_top = {}
            for query in query_labels:
                ranked = dist.loc[query, reference_labels].sort_values()
                tree_top[query] = ranked.index[0]
            cached_top = (
                cached_hits[cached_hits["rank"] == 1]
                .set_index("query_label")["reference_label"]
                .to_dict()
            )
            assert tree_top == cached_top, f"Computed closest hits do not match cached table: {tree_top} vs {cached_top}"

            closest = (
                cached_hits[cached_hits["rank"] == 1]
                .rename(columns={
                    "query_label": "query",
                    "reference_label": "closest_reference",
                    "fraction_different": "cached_hit_fraction_different",
                    "percent_identity_teaching_window": "percent_similarity",
                })
                [["query", "closest_reference", "cached_hit_fraction_different", "percent_similarity"]]
                .copy()
            )
            closest["tree_distance_to_closest"] = [
                float(dist.loc[row["query"], row["closest_reference"]])
                for _, row in closest.iterrows()
            ]
            display(closest.style.format({
                "cached_hit_fraction_different": "{:.4f}",
                "percent_similarity": "{:.2f}",
                "tree_distance_to_closest": "{:.4f}",
            }))
            print("Closest-reference labels match the cached hit table; identity percentages come from the cached direct-hit table.")
            """
        ),
        code(
            """
            #@title Distance matrix heatmap
            fig, ax = plt.subplots(figsize=(8, 6.5))
            sns.heatmap(
                dist,
                cmap="viridis",
                vmin=0,
                vmax=float(np.nanmax(dist.values)),
                square=True,
                linewidths=0.5,
                linecolor="white",
                cbar_kws={"label": "fraction of compared positions that differ"},
                ax=ax,
            )
            ax.set_title("Pairwise 16S marker distances", loc="left", fontsize=12)
            ax.set_xlabel("")
            ax.set_ylabel("")
            plt.xticks(rotation=35, ha="right")
            plt.yticks(rotation=0)
            plt.tight_layout()
            plt.show()
            """
        ),
        md(
            """
            ## 8. Build distance trees

            How does a distance matrix become a tree?

            UPGMA repeatedly joins the closest clusters. It is transparent, but it behaves like lineages changed at roughly similar rates.

            Neighbor joining also starts with pairwise distances, but it is less tied to that equal-rate assumption. Comparing both methods helps us see that a tree depends on the data and the method.
            """
        ),
        code(
            """
            #@title Build UPGMA and neighbor-joining trees
            def as_distance_matrix(dist_df):
                names = list(dist_df.index)
                lower = []
                for i, name in enumerate(names):
                    lower.append([float(dist_df.iloc[i, j]) for j in range(i + 1)])
                return DistanceMatrix(names, lower)

            dm = as_distance_matrix(dist)
            constructor = DistanceTreeConstructor()
            upgma_tree = constructor.upgma(dm)
            nj_tree = constructor.nj(dm)
            upgma_tree.rooted = True
            nj_tree.rooted = False

            print("UPGMA tree and neighbor-joining tree built from the same distance matrix.")
            """
        ),
        code(
            """
            #@title Plot the tree(s)
            def plot_tree(tree, title):
                fig, ax = plt.subplots(figsize=(10, 5.4))
                Phylo.draw(
                    tree,
                    axes=ax,
                    do_show=False,
                    show_confidence=False,
                    label_func=lambda clade: clade.name if clade.name else "",
                )
                ax.set_title(title, loc="left", fontsize=12)
                ax.set_xlabel("sequence-distance-derived branch length")
                ax.spines[["top", "right", "left"]].set_visible(False)
                ax.tick_params(axis="y", length=0)
                ax.grid(axis="x", color="#eeeeee", linewidth=0.8)
                plt.tight_layout()
                plt.show()

            if TREE_METHOD_TO_SHOW == "UPGMA":
                plot_tree(upgma_tree, "UPGMA distance tree")
            elif TREE_METHOD_TO_SHOW == "Neighbor joining":
                plot_tree(nj_tree, "Neighbor-joining distance tree")
            else:
                plot_tree(upgma_tree, "UPGMA distance tree")
                plot_tree(nj_tree, "Neighbor-joining distance tree")
            """
        ),
        md(
            """
            ## 9. Report the result carefully

            What exactly did we build?

            We built a distance-based gene tree from one 16S marker window. That can support a closest-reference statement. It does not prove exact species identity, and it is not a full species tree.
            """
        ),
        code(
            """
            #@title Create a student report sentence
            row = closest.set_index("query").loc[QUERY_TO_REPORT]
            report = (
                f"{QUERY_TO_REPORT} is closest to {row['closest_reference']} in this cached 16S marker comparison "
                f"({row['percent_similarity']:.2f}% similarity across the teaching window). "
                "Because this is one short 16S region, I can report a closest reference, not prove exact species identity or a complete species tree."
            )
            display(HTML(f'''
            <div style='max-width: 900px; border-left: 4px solid {OKABE_ITO["bluish_green"]}; padding: 10px 14px; background: #fafafa; font-size: 15px; line-height: 1.45;'>
              <b>Careful claim:</b><br>{report}
            </div>
            '''))
            """
        ),
        md(
            """
            ## 10. Optional: where IQ-TREE fits

            IQ-TREE is useful, but for a different teaching purpose.

            UPGMA and neighbor joining are distance methods. They help students see the bridge from sequence differences to tree geometry.

            IQ-TREE is a model-based maximum-likelihood tool. Use it after this lesson if you want students to compare a simple distance tree with a modern model-based tree. Do not make it the default class path for the first run.
            """
        ),
        code(
            """
            #@title Optional advanced note: MAFFT/IQ-TREE path is off by default
            RUN_ADVANCED_IQTREE_PATH = False #@param {type:"boolean"}
            if RUN_ADVANCED_IQTREE_PATH:
                print("Advanced path:")
                print("1. Align cleaned project reads and references with MAFFT.")
                print("2. Run IQ-TREE with model selection, for example: iqtree2 -s aligned.fasta -m MFP -B 1000 -T AUTO")
                print("3. Compare the maximum-likelihood tree with the UPGMA/NJ teaching trees.")
                print("Keep this optional so the main class run remains low-friction.")
            else:
                print("Advanced IQ-TREE path skipped. The class-safe UPGMA/NJ workflow is complete.")
            """
        ),
        md(
            """
            ## 11. Replace the teaching cache later

            How does this become your real soil microbiome project?

            Keep the notebook structure the same, but replace the teaching FASTA and metadata with your team's cleaned 16S reads and selected database references. The critical rule stays the same: prepare the cache before class, push it to GitHub, and let Colab load known files.
            """
        ),
        code(
            """
            #@title Project replacement schema
            project_schema = pd.DataFrame([
                {
                    "file": "project_16s_reads.fasta",
                    "required columns or fields": "FASTA id, DNA sequence",
                    "example": ">TeamA_ASV_001 sample=Rhizosphere_A",
                },
                {
                    "file": "project_16s_references.fasta",
                    "required columns or fields": "FASTA id, accession, source database",
                    "example": ">Bacillus_ref accession=NR_102783.2 source=NCBI",
                },
                {
                    "file": "project_16s_metadata.csv",
                    "required columns or fields": "label, role, accession, source_database, source_url, date_retrieved, taxonomy, note",
                    "example": "TeamA_ASV_001, query, blank, class sample, blank, 2026-05-25, unknown, cleaned ASV",
                },
                {
                    "file": "project_16s_abundance_table.csv",
                    "required columns or fields": "sample_id plus one column per ASV",
                    "example": "Rhizosphere_A, 128, 34, ...",
                },
            ])
            wrapped_table(project_schema, ["file", "required columns or fields", "example"])

            print("GitHub cache pattern:")
            print('CACHE_BASE_URL = "https://raw.githubusercontent.com/<org>/<repo>/main/soil_16s_class_cache"')
            print("Set USE_GITHUB_CACHE=True only after the cache folder is pushed.")
            """
        ),
        code(
            """
            #@title Export class outputs
            output_dir = Path("soil_microbiome_16s_outputs")
            output_dir.mkdir(exist_ok=True)
            dist.to_csv(output_dir / "soil_16s_distance_matrix.csv")
            closest.to_csv(output_dir / "soil_16s_closest_reference_report.csv", index=False)
            Phylo.write(upgma_tree, output_dir / "soil_16s_upgma_tree.newick", "newick")
            Phylo.write(nj_tree, output_dir / "soil_16s_neighbor_joining_tree.newick", "newick")
            metadata.to_csv(output_dir / "soil_16s_metadata_used.csv", index=False)
            print("Wrote outputs to:", output_dir.resolve())
            """
        ),
        md(
            """
            ## Final Think Prompts

            - Which query has the closest reference in the cached set?
            - Which references cluster near each other?
            - Does a high 16S similarity prove exact species identity?
            - What extra evidence would you want before making a stronger species claim?
            - How would the workflow change when you replace these cached teaching reads with your team's real project reads?
            """
        ),
    ]

    nb = {
        "cells": cells,
        "metadata": {
            "colab": {"name": NOTEBOOK_PATH.name, "provenance": []},
            "kernelspec": {"display_name": "Python 3", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    return nb


def validate_notebook(nb: dict) -> None:
    for index, cell in enumerate(nb["cells"]):
        if cell["cell_type"] == "code":
            ast.parse(cell["source"], filename=f"cell_{index}")


def fraction_different(seq_a: str, seq_b: str) -> float:
    n = min(len(seq_a), len(seq_b))
    compared = 0
    differences = 0
    for left, right in zip(seq_a[:n], seq_b[:n]):
        if left in "N-" or right in "N-":
            continue
        compared += 1
        if left != right:
            differences += 1
    return differences / compared if compared else math.nan


def smoke_test_core(cache_files: dict[str, str]) -> dict[str, object]:
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    from io import StringIO
    import xml.etree.ElementTree as ET

    refs = parse_fasta(cache_files["pilot_16s_references.fasta"])
    queries = parse_fasta(cache_files["pilot_16s_query_reads.fasta"])
    cached_hits = list(csv.DictReader(StringIO(cache_files["pilot_16s_cached_hits.csv"])))
    cached_blast_root = ET.fromstring(cache_files["pilot_16s_cached_blast.xml"])
    cached_xml_hits = cached_blast_root.findall(".//Hit")
    assert len(refs) == 5, f"expected 5 references, found {len(refs)}"
    assert len(queries) == 2, f"expected 2 queries, found {len(queries)}"
    assert len(cached_hits) == 10, f"expected 10 cached hit rows, found {len(cached_hits)}"
    assert len(cached_xml_hits) == 10, f"expected 10 XML hit rows, found {len(cached_xml_hits)}"

    windows: dict[str, str] = {}
    for label, (_header, seq) in refs.items():
        windows[label] = extract_marker(seq, 520)
    for label, (_header, seq) in queries.items():
        windows[label] = seq[:520]

    aligned_windows = star_align_to_reference(windows, "Bacillus_subtilis_168")
    labels = list(aligned_windows)
    distances: dict[tuple[str, str], float] = {}
    for left in labels:
        for right in labels:
            value = fraction_different(aligned_windows[left], aligned_windows[right])
            assert math.isfinite(value), f"non-finite distance for {left}, {right}"
            distances[(left, right)] = value

    def closest_reference(query_label: str) -> str:
        return min(refs, key=lambda ref_label: distances[(query_label, ref_label)])

    closest_a = closest_reference("Soil_ASV_A")
    closest_b = closest_reference("Soil_ASV_B")
    assert closest_a == "Bacillus_subtilis_168", closest_a
    assert closest_b == "Rhizobium_leguminosarum_IAM12609", closest_b
    cached_top = {
        row["query_label"]: row["reference_label"]
        for row in cached_hits
        if row["rank"] == "1"
    }
    assert cached_top == {"Soil_ASV_A": closest_a, "Soil_ASV_B": closest_b}, cached_top
    xml_top = {}
    for iteration in cached_blast_root.findall(".//Iteration"):
        query = iteration.findtext("Iteration_query-def")
        first_hit = iteration.find("./Iteration_hits/Hit[Hit_num='1']")
        xml_top[query] = first_hit.findtext("Hit_id")
    assert xml_top == cached_top, xml_top

    lower = []
    for i, left in enumerate(labels):
        lower.append([distances[(left, labels[j])] for j in range(i + 1)])
    dm = DistanceMatrix(labels, lower)
    constructor = DistanceTreeConstructor()
    upgma_tree = constructor.upgma(dm)
    nj_tree = constructor.nj(dm)
    assert len(upgma_tree.get_terminals()) == len(labels)
    assert len(nj_tree.get_terminals()) == len(labels)

    return {
        "references": len(refs),
        "queries": len(queries),
        "closest": {"Soil_ASV_A": closest_a, "Soil_ASV_B": closest_b},
        "cached_hit_rows": len(cached_hits),
        "cached_blast_xml_hits": len(cached_xml_hits),
        "terminal_count": len(labels),
    }


def main() -> None:
    CACHE_DIR.mkdir(exist_ok=True)
    cache_files = make_cache_files()
    for name, content in cache_files.items():
        (CACHE_DIR / name).write_text(content, encoding="utf-8", newline="\n")

    nb = make_notebook(cache_files)
    validate_notebook(nb)
    smoke = smoke_test_core(cache_files)
    (CACHE_DIR / "cache_validation_report.json").write_text(
        json.dumps(
            {
                "status": "passed",
                "checks": smoke,
                "note": "Generated by build_soil_microbiome_16s_class_safe_pilot.py",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    NOTEBOOK_PATH.write_text(json.dumps(nb, indent=2), encoding="utf-8", newline="\n")

    print(f"Wrote {NOTEBOOK_PATH}")
    print(f"Wrote cache files under {CACHE_DIR}")
    print(f"Validated {sum(1 for c in nb['cells'] if c['cell_type'] == 'code')} code cells with ast.parse")
    print(f"Smoke-tested core workflow: {json.dumps(smoke, sort_keys=True)}")


if __name__ == "__main__":
    main()
