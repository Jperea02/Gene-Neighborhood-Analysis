#!/usr/bin/env python3
"""
Yarrowia engineered-strain genomic context pipeline
===================================================

Single-file, step-by-step workflow for W29/PO1f-derived engineered Yarrowia strains.

What this pipeline does
-----------------------
STEP 1  - Build motif FASTA from sequences embedded in this script
STEP 2  - Load strain contigs and CLIB122 reference FASTA
STEP 3  - Assign contigs to CLIB122 chromosomes / mtDNA by k-mer similarity
STEP 4  - Detect engineered motifs on contigs by exact sequence matching
STEP 5  - Group motif hits into candidate integration loci
STEP 6A - Build contig-to-reference structural alignment blocks with minimap2 (PAF)
STEP 6B - Infer insertion intervals using paired flank alignments constrained by structural blocks
STEP 7  - Parse CLIB122 GFF3 into a clean annotation table
STEP 8  - Build mapped essential-gene table (CLIB89 -> CLIB122)
STEP 9  - Extract instability-associated features from CLIB122 annotation
STEP 10 - Perform confidence-aware neighborhood analysis:
          genes within insertion region, nearest flanking genes,
          essential genes in window, and nearby instability features
STEP 11 - Write summary tables

Important note
--------------
Your strain is PO1f-derived, so CLIB89 essentiality is biologically closer to your strain
background. Structural placement of insertion loci is inferred in CLIB122 reference space
using contig alignment blocks and flank-based breakpoint mapping. Therefore:

- high-confidence loci can support exact insertion-site interpretation
- moderate-confidence loci should be interpreted as localized insertion regions
- distances to essential genes should be interpreted as distances to the nearest
  CLIB122 ortholog of a CLIB89 essential gene

Current motif detection
-----------------------
This version uses exact motif matching on both strands. That works well for polished
contigs and known engineered sequences. If you later want mismatch-tolerant matching,
you can upgrade this to BLAST or minimap2-based motif alignment.

Dependencies
------------
Required:
    pandas
    biopython
    openpyxl   (if your tables are .xlsx)
    minimap2

Optional:
    samtools   (not required by this script, but useful if you want BAM files later)

Install Python deps:
    pip install pandas biopython openpyxl

Install minimap2 on Mac:
    brew install minimap2

Install minimap2 with conda:
    conda install -c bioconda minimap2
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
from pathlib import Path
from textwrap import wrap
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO


# =============================================================================
# USER CONFIGURATION
# =============================================================================

CONFIG = {
    # ----------------------------
    # Project / strain
    # ----------------------------
    "strain_name": "YL-W3-FAT-FAA",
    "contigs_fasta": "YL-W3-FAT-FAA_reference.fasta",

    # ----------------------------
    # Reference files
    # ----------------------------
    "clib122_genome_fasta": "../reference_data/CLIB122_genome.fna",
    "clib122_gff3": "../reference_data/CLIB122_annotation.gff",
    "clib89_essential_table": "../reference_data/CLIB89_essential_genes.xlsx",
    "clib89_clib122_mapping_table": "../reference_data/CLIB89_CLIB122_mapping.xlsx",

    # ----------------------------
    # Output folder
    # ----------------------------
    "output_dir": "results_YL-W3-FAT-FAA",

    # ----------------------------
    # Core settings
    # ----------------------------
    "integration_grouping_distance_bp": 10000,
    "gene_neighborhood_n": 5,
    "essential_neighborhood_n": 5,
    "instability_neighborhood_n": 5,

    # ----------------------------
    # K-mer assignment settings
    # ----------------------------
    "kmer_size": 17,

    # ----------------------------
    # minimap2 upgrade
    # ----------------------------
    "use_minimap2_projection": True,
    "minimap2_preset": "asm5",   # good for polished assembly-to-reference mapping
}

CHR_MAP = {
    "NC_006067.1": "ChrA",
    "NC_006068.1": "ChrB",
    "NC_006069.1": "ChrC",
    "NC_006070.1": "ChrD",
    "NC_006071.1": "ChrE",
    "NC_006072.1": "ChrF",
    "NC_002659.1": "mtDNA",
}

# Put your engineered motif / construct sequences here.
# Replace the placeholder sequences with your real sequences.
#
# You can include:
# - CDSs (MhFAR, GFP, URA3)
# - promoters
# - terminators
# - linkers
# - fusion junctions
#
# REPLACE with real motif sequences
MOTIFS = {
    "MhFAR": "gccatccagcaggtccaccacgccgacacctcctcctccaaggtcctgggacagctgcgaggaaagcgagtcctgatcaccggcaccaccggcttcctgggcaaggtcgtcctggagcggctgattcgagctgtccctgacatcggagctatctatctgctgatccggggcaacaagcgacaccctgacgctcgatctcgattcctggaggagatcgccacctcctccgtcttcgaccggctgcgagaggctgactctgagggattcgacgccttcctggaggagcggatccactgcgtcaccggcgaggtcaccgaggctggattcggcatcggacaggaggactatcggaagctggccaccgagctggacgctgtcatcaactccgccgcctccgtcaacttccgggaggagctggacaaggccctggccatcaacaccctgtgcctgcggaacatcgccggcatggtcgacctgaaccccaagctggccgtcctgcaggtctccacctgctatgtcaacggcatgaactccggccaggtcaccgagtccgtcatcaagcctgctggagaggctgttcctcgttctcctgacggcttctatgagatcgaggagctggtccggctgctgcaggacaagatcgaggacgtccaggcccggtattccggcaaggtcctggagcggaagctggtcgacctgggcatccgggaggccaaccggtatggctggtccgacacctataccttcaccaagtggctgggcgagcagctgctgatgaaggccctgaacggacgaaccctgaccatcctgcggccctccatcatcgagtccgctctggaggagcctgctcctggctggatcgagggcgtcaaggtcgccgacgccatcatcctggcctatgcccgggagaaggtcaccctgttccccggcaagcggtctggaatcatcgacgtcatccccgtcgacctggtcgccaactccatcatcctgtccctggctgaggctctgggagagcctggacgacgacgaatctatcagtgctgctccggaggaggcaaccctatctccctgggcgagttcatcgaccacctgatggccgagtccaaggccaactatgccgcctatgaccacctgttctatcggcagccctccaagcccttcctggccgtcaaccgagctctgttcgacctggtcatctccggcgtccgactgcctctgtccctgaccgaccgggtcctgaagctgctgggcaactcccgggacctgaagatgctgcggaacctggacaccacccagtccctggccaccatcttcggcttctataccgctcctgactatatcttccggaacgacgagctgatggccctggccaaccggatgggcgaggtcgacaagggcctgttccccgtcgacgcccgactgatcgactgggagctgtatctgcggaagatccacctggccggcctgaaccggtatgccctgaaggagcggaaggtctattccctgaagaccgcccggcagcggaagaaggccgcctga",
    "GFP": "atggtgagcaagcagatcctgaagaacaccggcctgcaggagatcatgagcttcaaggtgaacctggagggcgtggtgaacaaccacgtgttcaccatggagggctgcggcaagggcaacatcctgttcggcaaccagctggtgcagatccgcgtgaccaagggcgcccccctgcccttcgccttcgacatcctgagccccgccttccagtacggcaaccgcaccttcaccaagtaccccgaggacatcagcgacttcttcatccagagcttccccgccggcttcgtgtacgagcgcaccctgcgctacgaggacggcggcctggtggagatccgcagcgacatcaacctgatcgaggagatgttcgtgtaccgcgtggagtacaagggccgcaacttccccaacgacggccccgtgatgaagaagaccatcaccggcctgcagcccagcttcgaggtggtgtacatgaacgacggcgtgctggtgggccaggtgatcctggtgtaccgcctgaacagcggcaagttctacagctgccacatgcgcaccctgatgaagagcaagggcgtggtgaaggacttccccgagtaccacttcatccagcaccgcctggagaagacctacgtggaggacggcggcttcgtggagcagcacgagaccgccatcgcccagctgaccagcctgggcaagcccctgggcagcctgcacgagtgggtgtaa",
    "URA3": "atgccctcctacgaagctcgagctaacgtccacaagtccgcctttgccgctcgagtgctcaagctcgtggcagccaagaaaaccaacctgtgtgcttctctggatgttaccaccaccaaggagctcattgagcttgccgataaggtcggaccttatgtgtgcatgatcaagacccatatcgacatcattgacgacttcacctacgccggcactgtgctccccctcaaggaacttgctcttaagcacggtttcttcctgttcgaggacagaaagttcgcagatattggcaacactgtcaagcaccagtacaagaacggtgtctaccgaatcgccgagtggtccgatatcaccaacgcccacggtgtacccggaaccggaatcattgctggcctgcgagctggtgccgaggaaactgtctctgaacagaagaaggaggacgtctctgactacgagaactcccagtacaaggagttcctggtcccctctcccaacgagaagctggccagaggtctgctcatgctggccgagctgtcttgcaagggctctctggccactggcgagtactccaagcagaccattgagcttgcccgatccgaccccgagtttgtggttggcttcattgcccagaaccgacctaagggcgactctgaggactggcttattctgacccccggggtgggtcttgacgacaagggagacgctctcggacagcagtaccgaactgttgaggatgtcatgtctaccggaacggatatcataattgtcggccgaggtctgtacggccagaaccgagatcctattgaggaggccaagcgataccagaaggctggctgggaggcttaccagaagattaactgttag",
    "Hyg": "tcactccttggctctgggtctggtggagggtcttctgtttccagagtcagccaggacctcgacgcatccgtcggtccagacagcagcggatcttctggcgatctgggttcttccgacggttccagctccggatctgacgatggcgtcgcatcttccctgagcccaagcagcgtcgtcgaagttgccgtcgaccagggactgatacagctggtccaggccgattctcagcatataggctctcagtctaggagatccagccagctcggggtgtcttctctcgaaatatctggtctgctgctccatgcaagccagccagggtctccagaagaagatgttggcgacctcatactgggagtcgccgaacatggcctcggaccagtcgatgacggcggtgattctgccgttgtcggtcaggacgttgttggagccgaagtcggcgtggaccaggtgtctgacctcaggacagtcctcggcccacagcatcagctcgtccagagcctgagcgacagaagcggagacggtgtcgtccatgacggtctgccagtgatagacatgaggatcggcgatggcgcagatgaagtctctccaggtggtatactggccgatgccttgaggtccgaaaggtccgaagccggaggtctgagacaggtcagcagcggcgatagcgtccattgcctctgcgacaggctgcaggacagcaggcagctcggtctcaggcaggtcctgcagagtgactccctgggctcttctggagatgcaataggtcagggactcggagaactcgccgatgtccaggacctcagggataggcagagcagcagaggcgaagtgccgatagacgtaccggtccttatagaagccgtcggcgcaggagttgactctcaggacatatcctctgccgccgacgtcgaaggagaaggctctagactcctctccctcggacagctgcatcaggtcggagacggagtcgaacttctcgatcaggaacttctcgacagaagtggcggtcagctcgggcttcttcat",
    "WS": "cgacccctccaccctattgatttcattttcctctctctggagaagcgacagcagcctatgcacgttggcggactcttcctctttcagatccccgacaacgctcctgataccttcattcaggacctggttaacgacatccgaatttcgaagtccattcccgtccccccttttaacaacaagctgaacggcctcttctgggacgaggatgaggagtttgacctggatcaccatttccgacacatcgccctcccccatcctggtcgaatccgagagctgctcatctacatttcgcaggagcactccaccctgctcgaccgagctaagcccctgtggacttgcaacatcattgagggcatcgagggtaaccgattcgccatgtactttaagattcaccatgctatggtggacggagttgccggcatgcgactgattgagaagtctctctcgcatgacgtcaccgagaagtccatcgttcccccttggtgtgtcgagggaaagcgagctaagcgactgcgagagcccaagactggcaagatcaagaagattatgtctggtatcaagtcgcagctgcaggccacccctactgtgattcaggagctgtcccagaccgtttttaaggacatcggccgaaaccccgatcacgtctcctctttccaggccccttgctctatcctgaaccagcgagtgtcgtcctctcgacgattcgccgctcagtcttttgacctcgatcgattccgaaacatcgctaagtcgctgaacgtcactattaacgacgtggtgctcgccgtgtgttctggtgccctgcgagcttacctcatgtcccacaactctctgccctcgaagcctctcatcgctatggtccccgcttctattcgaaacgacgattccgacgtttctaaccgaatcaccatgattctggccaacctcgctactcataaggacgatcctctgcagcgactggagatcattcgacgatccgtgcagaactctaagcagcgattcaagcgaatgacctcggaccagattctgaactactccgctgttgtctacggtcccgccggactcaacatcatttccggaatgatgcctaagcgacaggctttcaacctggtcatctctaacgtgcctggaccccgagagcctctctactggaacggtgctaagctggacgctctctaccctgcctccatcgttctggatggtcaggccctcaacattaccatgacttcttacctggacaagctggaggtcggactgatcgcttgccgaaacgccctcccccgaatgcagaacctgctcacccacctggaagaagagattcagctcttcgaaggtgttattgccaagcaggaagacatcaagactgccaactaa",
    "FAT1": "ATGAAAACGATATTGAAAATAACAAAATCCGAGAATCAAAACGCACTTTTCAAAAATCCCATATCTCCTCCACATCCCCCACAAACCCGAACCCCTAGTCTCAAAATTAAGGTTCAACCCCAGATACCCCATTTTTTTCATGCTGGCCCTTATATTAACCGTGGCTGCCCTTTTTTAAGTCCCCTTCTCCACTATCACCTTGTCGAAATCCCAACCACAATGACAGCTGGACTAGTTGCTGCCGCCGCCATTGGAGCCGCCTACCTGGAGGCCAAGACGCTCATCTCCGAAGATGCCTACATGATCCGAGGCGCCATGACCAACGGTCTGGATTTTTTCTACAACGCCTGGAAGGGCCGAGTGCAATACTGGTACGCGTTTGAAGACGCGGTGAAGAAGTACCCCAACAACCCGGCAATCGTCTACCCCAAGCCCATCGAGGGCAAGAAGCCCAGCGGCGACTCGTACGACGACCTGTTTGACGTGGAGACCTTCACCTACCAGCAGCTGTACGATGAGGTGCTGAAAATGTCACACCTGCTGCGAAACAAGTACGGAGTGACCGCCAACGACACCATTGCCCTCAACGCCATGAACTCGCCCCTGTTCATTATCGTGTGGTTCGCCATCTGGAACCTGGGCGCCACCCCGGCCTTCATTAACTACAACCTGGCCGACAAGTCGCTGCTGCACTGCCTCAAGGTCGGACATGCCTCCATCATGTTCGTCGACACCGAGGTCGAGGGCAACGTGCGACCTTCGCTGGCCGAGATCAAGTCCGAGGCCAAGTGCGACACCGTGTTCATGGACGACGACTTCCTCGCCGCCTACGCCGCCTCCCCGGCCTACAGAGCCCCCGACTACGAGCGACATCCCGAGCAGAAGGACTATGACACTGCTGTGCTCATCTACACGTCCGGAACCACCGGTCTGCCCAAGCCCGCCATCATGTCGTGGAAGAAGGCCAAGCTCATGTCGTCGCTCTACGGCCACTCCATCCGGCTCAAGAACAACGGTGTCGTCTACTCCGCCATGCCTCTCTACCACTCCACCGCCGCCATTCTCGGCTGTCTGCCCTGCCTCAACCGAGGAGCCGCTTACGCCCCTGGCCGAAAGTTCTCTACCACCACCTTCTGGACCCAGGCCAAGCTCACCAACGCCACTCACATCCAGTACGTCGGAGAGACCTGCCGATACCTGATCAACGCCCCTCCCTCTCCCGACGAGAAGAGTCACCAGATCAAGGTGGCGTTCGGCAACGGAATGCGACGAGATATTTGGGTCAAGTTCAAGGAACGATTCAACATCCCCGCCATCGGCGAGTTTTACGCCGCCACCGAGGGTCCTCTCGGAACCAACAACTTCCAGCAGGGTGAGATTGGTATCGGTGCCATGGGCCGATACGGTAAGCTTCTGGCTGCCATCCTGGCCACCCGACAGACCATCGTTCCTGTGGATCCCGAAGACGAGACTGAGCTGTGGCGAGACCCCGAGACCGGCTTCTGCCGAGTCGCCCAGTCCGACGAGCCCGGTGAGTTCATCCAGAAGATCCCCAACCCCGAGAAGGTGCACGAGACCTTCCAGGGATACCTTGGTAACGATAAGGCCACCAACTCCAAGATTATGCGAGACGTGTTCAAGAAGGGCGATGCCTACTACCGAACCGGTGATCTGGTGCGACTCAACGACGAGCAGTGCTACTACTTTGTCGACCGTCTTGGAGACACCTTCCGATGGAAGTCCGAGAATGTGTCCACCTCCGAGGTGGAGGAGCACGTTGGCGCTTCCGACCCCAACATTGAGCAGGTTGTCTGCGTGGGTGTCAAGGTGCCCGAGCACGAGGGCCGAGCTGGTTTCGCCGTGGTGAAGCTCAAGGACGCTTCTGTCAAGCCCAACCTTGACCAGATTGCCGAGTACTCGCTCAAGCAGCTGCCCAAGTACGCCGTGCCTCTGTTCATCAAGTTTGTCGACGAGATTGAGCGAACCGGCAACAACAAGGTGCAGAAGGTCAAGTACAAGAACCAGAAGATGCCCCATGAGGAAGGCGAGAGCCCCATTTACTGGCTCAAGGGCAACAAGTACGTGGAGCTCGATGCTGGAGACTGGGCTTCTCTGGGATCCGGAAAGATTAAGCTGTAA",
    "FAA1": "ATGGTCGGATACACAATTTCCTCAAAGCCCGTGTCGGTGGAGGTCGGCCCCGCCAAGCCTGGCGAGACTGCCCCCCGACGAAACGTCATTGCCAAGGACGCCCCTGTCGTCTTCCCCGACAACGACTCGTCCCTGACCACCGTCTACAAGCTGTTCAAAAAGTACGCCGAGATCAACAGCGAGCGAAAGGCCATGGGATGGCGAGACACCATCGACATCCACGTGGAGACCAAACAGGTGACCAAGGTCGTGGACGGAGTGGAGAAGAAGGTGCCCAAGGAATGGAAGTACTTTGAGATGGGCCCTTACAAGTGGCTCTCATACAAGGAGGCCCTTAAGCTGGTCCATGATTATGGAGCTGGTCTTCGACACCTCGGAATCAAGCCCAAGGAGAAGATGCACATTTACGCCCAGACCTCCCACCGATGGATGCTCTCTGGCCTGGCTTCTCTGTCTCAGGGTATTCCCATTGTCACTGCCTACGACACTCTTGGAGAGGAGGGTCTCACTCGATCTCTCCAGGAGACCAACTCGGTCATCATGTTTACCGACAAGGCTCTGCTGAGCTCTCTCAAGGTCTCTCTCAAGAAGGGCACCGATCTGCGAATCATCATCTACGGAGGTGATCTGACCCCCGACGACAAGAAGGCCGGAAACACGGAGATTGACGCCATCAAGGAGATTGTTCCAGATATGAAGATCTACACCATGGACGAGGTTGTCGCTCTCGGCCGAGAACACCCCCACCCCGTGGAGGAGGTCGACTATGAGGACCTGGCCTTCATCATGTACACCTCTGGTTCTACCGGTGTCCCCAAGGGTGTGGTTCTGCAGCACAAGCAGATCCTCGCCTCTGTGGCCGGTGTCACCAAGATCATTGACCGATCTATCGTCGGCAACACAGACCGGCTTCTCAACTTCCTGCCCCTCGCACACATTTTCGAGTTTGTGTTCGAGATGGTCACCTTCTGGTGCGGTGCTTCTCTGGGTTACGGAACCGTCAAGACCATTTCCGATCTGTCCATGAAGAACTGTAAGGGAGACATTCGAGAGCTCAAGCCCACCATCATGGTCGGCGTTCCCGCTGTCTGGGAACCTATGCGAAAGGGTATTCTTGGCAAGATCAAGGAGCTGTCTCCTCTGATGCAGCGGGTCTTCTGGGCCTCATTTGCCGCCAAGCAGCGTCTCGACGAGAACGGACTCCCTGGTGGATCTATCCTCGACTCGCTCATTTTCAAGAAGGTCAAGGACGCCACTGGAGGCTGTCTCCGATACGTGTGTAACGGAGGTGCTCCAGTATCTGTCGACACCCAGAAGTTCATCACCACTCTCATCTGTCCCATGCTGATTGGATGCGGTCTGACCGAGACTACAGCCAACACCACCATCATGTCGCCTAAATCGTACGCCTTTGGCACCATTGGTGAGCCCACCGCCGCCGTGACCCTCAAGCTCATTGACGTGCCTGAAGCCGGCTACTTCGCCGAGAACAACCAGGGAGAGCTGTGCATCAAGGGCAACGTCGTGATGAAGGAGTACTACAAGAACGAGGAGGAGACCAAGAAGGCGTTCTCCGACGATGGCTACTTCCTCACCGGTGATATTGCCGAGTGGACCGCCAATGGCCAGCTCAGAATCATTGACCGACGAAAGAACCTCGTCAAGACCCAGAACGGAGAGTACATTGCTCTGGAGAAGCTCGAGACACAGTACCGATCGTCGTCGTACGTGGCCAACCTGTGTGTGTACGCCGACCAGAACCGAGTCAAGCCCATTGCTCTGGTCATTCCTAACGAGGGCCCCACCAAGAAGCTTGCCCAGAGCTTGGGCGTCGATTCTGACGACTGGGACGCCGTCTGCTCCAACAAAAAGGTGGTCAAGGCTGTGCTCAAGGACATGCTCGATACCGGCCGATCTCTGGGTCTGTCCGGCATTGAGCTGCTGCAAGGCATTGTGTTGCTGCCTGGCGAGTGGACTCCTCAGAACAGCTACCTGACTGCTGCCCAGAAGCTCAACCGAAAGAAGATTGTGGATGATAACAAGAAGGAAATTGATGAGTGCTACGAGCAGTCTTAG",
}


# =============================================================================
# LOGGING / SMALL UTILITIES
# =============================================================================

def step(msg: str) -> None:
    print("\n" + "=" * 88)
    print(msg)
    print("=" * 88)

def info(msg: str) -> None:
    print(f"[INFO] {msg}")

def warn(msg: str) -> None:
    print(f"[WARN] {msg}")

def fail(msg: str) -> None:
    raise RuntimeError(msg)

def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p

def clean_seq(seq: str) -> str:
    return "".join(str(seq).upper().split())

def validate_dna(seq: str, name: str = "sequence") -> None:
    bad = sorted(set(seq) - set("ACGTN"))
    if bad:
        fail(f"{name} contains invalid DNA characters: {bad}")

def rc(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]

def fasta_write(records: Dict[str, str], out_fasta: str | Path) -> None:
    out_fasta = Path(out_fasta)
    with out_fasta.open("w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for line in wrap(seq, 80):
                fh.write(line + "\n")

def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if suffix in {".csv", ".tsv", ".txt"}:
        sep = "\t" if suffix == ".tsv" else ","
        return pd.read_csv(path, sep=sep, low_memory=False)
    fail(f"Unsupported table format: {path}")

def normalize_yali0(x: str) -> str:
    x = str(x).strip()
    if re.match(r"^YALI0_[A-Z]\d+[a-z]?$", x):
        return x
    m = re.match(r"^(YALI0)([A-Z]\d+[a-z]?)$", x)
    if m:
        return f"{m.group(1)}_{m.group(2)}"
    return x

def load_fasta_as_dict(path: str | Path) -> Dict[str, str]:
    path = Path(path)
    if not path.exists():
        fail(f"Missing FASTA file: {path}")
    records = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        records[rec.id] = clean_seq(str(rec.seq))
    if not records:
        fail(f"No sequences found in FASTA: {path}")
    return records

def extract_kmers(seq: str, k: int) -> set[str]:
    seq = clean_seq(seq)
    if len(seq) < k:
        return set()
    return {seq[i:i+k] for i in range(len(seq) - k + 1) if "N" not in seq[i:i+k]}

def jaccard(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)

def interval_distance(a_start: int, a_stop: int, b_start: int, b_stop: int) -> int:
    if a_stop >= b_start and a_start <= b_stop:
        return 0
    if a_stop < b_start:
        return b_start - a_stop
    return a_start - b_stop

def detect_headerless_ids(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    for c in df.columns:
        lc = c.lower()
        for t in candidates:
            if t.lower() == lc:
                return c
    return None


# =============================================================================
# STEP 1 - BUILD MOTIF FASTA
# =============================================================================

def build_motif_fasta(out_dir: Path) -> Tuple[Path, Dict[str, str]]:
    step("STEP 1 — Build motif FASTA from sequences embedded in this Python file")
    cleaned = {}
    for name, seq in MOTIFS.items():
        seq2 = clean_seq(seq)
        validate_dna(seq2, name)
        cleaned[name] = seq2

    out_fasta = out_dir / "motifs.fasta"
    fasta_write(cleaned, out_fasta)

    info(f"Saved motif FASTA: {out_fasta}")
    for name, seq in cleaned.items():
        info(f"  {name}: {len(seq)} bp")

    return out_fasta, cleaned


# =============================================================================
# STEP 2 - LOAD FASTAS
# =============================================================================

def load_inputs() -> Tuple[Dict[str, str], Dict[str, str]]:
    step("STEP 2 — Load strain contigs and CLIB122 reference FASTA")
    contigs = load_fasta_as_dict(CONFIG["contigs_fasta"])
    refs = load_fasta_as_dict(CONFIG["clib122_genome_fasta"])
    info(f"Loaded {len(contigs)} contig(s) from {CONFIG['contigs_fasta']}")
    info(f"Loaded {len(refs)} reference record(s) from {CONFIG['clib122_genome_fasta']}")
    return contigs, refs


# =============================================================================
# STEP 3 - CONTIG -> REFERENCE ASSIGNMENT BY K-MERS
# =============================================================================

def assign_contigs_to_reference(contigs: Dict[str, str], refs: Dict[str, str], out_dir: Path) -> pd.DataFrame:
    step("STEP 3 — Assign contigs to CLIB122 chromosomes / mtDNA by k-mer similarity")

    k = int(CONFIG["kmer_size"])
    ref_kmers = {name: extract_kmers(seq, k) for name, seq in refs.items()}
    rows = []

    for contig_name, contig_seq in contigs.items():
        ck = extract_kmers(contig_seq, k)
        best_name = None
        best_score = -1.0

        for ref_name, ref_seq in refs.items():
            score = jaccard(ck, ref_kmers[ref_name])
            if score > best_score:
                best_score = score
                best_name = ref_name

        if best_score >= 0.5:
            assigned_ref = best_name
            assigned_label = CHR_MAP.get(best_name, best_name)
            assignment_confidence = "high"
        elif best_score >= 0.05:
            assigned_ref = best_name
            assigned_label = CHR_MAP.get(best_name, best_name)
            assignment_confidence = "moderate"
        else:
            assigned_ref = "UNRESOLVED"
            assigned_label = "UNRESOLVED"
            assignment_confidence = "low"

        rows.append({
            "contig": contig_name,
            "assigned_reference": assigned_ref,
            "assigned_reference_label": assigned_label,
            "kmer_jaccard": best_score,
            "assignment_confidence": assignment_confidence,
            "contig_length": len(contig_seq),
        })

    df = pd.DataFrame(rows).sort_values(["assigned_reference", "contig"])
    out = out_dir / "contig_assignment.csv"
    df.to_csv(out, index=False)
    info(f"Saved: {out}")
    print(df.to_string(index=False))
    return df

# =============================================================================
# STEP 4 - EXACT MOTIF DETECTION
# =============================================================================

def find_all_exact_hits(sequence: str, query: str) -> List[int]:
    hits = []
    start = 0
    while True:
        idx = sequence.find(query, start)
        if idx == -1:
            break
        hits.append(idx + 1)  # 1-based
        start = idx + 1
    return hits

def detect_motifs_on_contigs(contigs: Dict[str, str], motifs: Dict[str, str], out_dir: Path) -> pd.DataFrame:
    step("STEP 4 — Detect engineered motifs on contigs by exact sequence matching")

    rows = []
    for contig_name, contig_seq in contigs.items():
        for motif_name, motif_seq in motifs.items():
            motif_rc = rc(motif_seq)

            for pos in find_all_exact_hits(contig_seq, motif_seq):
                rows.append({
                    "Name": motif_name,
                    "Sequence Name": contig_name,
                    "Minimum": pos,
                    "Maximum": pos + len(motif_seq) - 1,
                    "Length": len(motif_seq),
                    "Direction": "forward",
                })

            for pos in find_all_exact_hits(contig_seq, motif_rc):
                rows.append({
                    "Name": motif_name,
                    "Sequence Name": contig_name,
                    "Minimum": pos,
                    "Maximum": pos + len(motif_seq) - 1,
                    "Length": len(motif_seq),
                    "Direction": "reverse",
                })

    hits = pd.DataFrame(rows)
    if hits.empty:
        warn("No motif hits were found. Check your motif sequences or contig assembly.")
    else:
        hits = hits.sort_values(["Sequence Name", "Minimum", "Maximum", "Name"])

    out = out_dir / "motif_hits.csv"
    hits.to_csv(out, index=False)
    info(f"Saved: {out}")
    if not hits.empty:
        print(hits.to_string(index=False))
    return hits


# =============================================================================
# STEP 5 - GROUP MOTIF HITS INTO LOCI
# =============================================================================

def group_hits_into_loci(hits: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    step("STEP 5 — Group motif hits into candidate integration loci")

    if hits.empty:
        empty = pd.DataFrame()
        empty.to_csv(out_dir / "integration_loci.csv", index=False)
        return empty

    max_gap = int(CONFIG["integration_grouping_distance_bp"])
    loci = []
    locus_counter = 1

    for contig, sub in hits.groupby("Sequence Name", sort=False):
        sub = sub.sort_values("Minimum").reset_index(drop=True)

        current_start = int(sub.loc[0, "Minimum"])
        current_stop = int(sub.loc[0, "Maximum"])
        current_rows = [sub.loc[0]]

        for i in range(1, len(sub)):
            row = sub.loc[i]
            row_start = int(row["Minimum"])
            row_stop = int(row["Maximum"])

            if row_start - current_stop <= max_gap:
                current_stop = max(current_stop, row_stop)
                current_rows.append(row)
            else:
                motif_names = sorted({r["Name"] for r in current_rows})
                strands = sorted({r["Direction"] for r in current_rows})
                loci.append({
                    "locus_id": f"L{locus_counter}",
                    "contig": contig,
                    "locus_start": current_start,
                    "locus_stop": current_stop,
                    "motifs_present": ",".join(motif_names),
                    "strand_set": ",".join(strands),
                    "motif_count": len(current_rows),
                })
                locus_counter += 1
                current_start = row_start
                current_stop = row_stop
                current_rows = [row]

        motif_names = sorted({r["Name"] for r in current_rows})
        strands = sorted({r["Direction"] for r in current_rows})
        loci.append({
            "locus_id": f"L{locus_counter}",
            "contig": contig,
            "locus_start": current_start,
            "locus_stop": current_stop,
            "motifs_present": ",".join(motif_names),
            "strand_set": ",".join(strands),
            "motif_count": len(current_rows),
        })
        locus_counter += 1

    loci_df = pd.DataFrame(loci)
    out = out_dir / "integration_loci.csv"
    loci_df.to_csv(out, index=False)
    info(f"Saved: {out}")
    print(loci_df.to_string(index=False))
    return loci_df


# =============================================================================
# STEP 6A - BUILD ALIGNMENT BLOCKS (PAF)
# =============================================================================

def run_minimap2_paf(query_fasta: str | Path, ref_fasta: str | Path, out_paf: str | Path) -> None:
    """
    Run minimap2 and output PAF format (used for alignment blocks + flank mapping)
    """
    if shutil.which("minimap2") is None:
        fail("minimap2 is not installed or not on PATH.")

    cmd = [
        "minimap2",
        "-x", "asm5",
        str(ref_fasta),
        str(query_fasta),
    ]

    info(f"Running minimap2 (PAF): {' '.join(cmd)}")

    with Path(out_paf).open("w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, check=True)

def parse_paf_blocks(paf_file: str | Path) -> pd.DataFrame:
    """
    Parse minimap2 PAF into an alignment block table.
    PAF columns:
      0 qname
      1 qlen
      2 qstart
      3 qend
      4 strand
      5 tname
      6 tlen
      7 tstart
      8 tend
      9 nmatch
     10 alnlen
     11 mapq
    """
    rows = []

    with Path(paf_file).open() as fh:
        for line in fh:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2]) + 1
            qend = int(parts[3])
            strand = parts[4]
            tname = parts[5]
            tlen = int(parts[6])
            tstart = int(parts[7]) + 1
            tend = int(parts[8])
            nmatch = int(parts[9])
            alnlen = int(parts[10])
            mapq = int(parts[11])

            pid = (nmatch / alnlen * 100) if alnlen > 0 else 0.0

            rows.append({
                "contig": qname,
                "contig_length": qlen,
                "contig_start": qstart,
                "contig_stop": qend,
                "ref_chr": tname,
                "ref_length": tlen,
                "ref_start": tstart,
                "ref_stop": tend,
                "strand": strand,
                "matches": nmatch,
                "aln_block_len": alnlen,
                "percent_identity": pid,
                "mapq": mapq,
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df["chr_label"] = df["ref_chr"].map(CHR_MAP).fillna(df["ref_chr"])
    return df

def filter_alignment_blocks(blocks_df: pd.DataFrame, min_block_len: int = 5000, min_pid: float = 90.0) -> pd.DataFrame:
    if blocks_df.empty:
        return blocks_df

    out = blocks_df[
        (blocks_df["aln_block_len"] >= min_block_len) &
        (blocks_df["percent_identity"] >= min_pid)
    ].copy()

    out = out.sort_values(
        ["contig", "mapq", "percent_identity", "aln_block_len"],
        ascending=[True, False, False, False]
    ).reset_index(drop=True)

    return out

def assign_loci_to_blocks(loci_df: pd.DataFrame, blocks_df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign each cassette locus to the structural alignment block that overlaps
    its contig interval best. If no block overlaps, choose the nearest block in
    contig space.
    """
    results = []

    if loci_df.empty or blocks_df.empty:
        return pd.DataFrame()

    for _, locus in loci_df.iterrows():
        contig = locus["contig"]
        locus_start = int(locus["locus_start"])
        locus_stop = int(locus["locus_stop"])

        sub = blocks_df[blocks_df["contig"] == contig].copy()
        if sub.empty:
            continue

        overlaps = []
        for _, blk in sub.iterrows():
            ov_start = max(locus_start, int(blk["contig_start"]))
            ov_stop = min(locus_stop, int(blk["contig_stop"]))
            ov = max(0, ov_stop - ov_start + 1)
            overlaps.append(ov)

        sub["locus_overlap_bp"] = overlaps

        if (sub["locus_overlap_bp"] > 0).any():
            best = sub.sort_values(
                ["locus_overlap_bp", "mapq", "percent_identity", "aln_block_len"],
                ascending=[False, False, False, False]
            ).iloc[0]
        else:
            dists = []
            for _, blk in sub.iterrows():
                d = interval_distance(
                    locus_start, locus_stop,
                    int(blk["contig_start"]), int(blk["contig_stop"])
                )
                dists.append(d)
            sub["block_distance_bp"] = dists
            best = sub.sort_values(
                ["block_distance_bp", "mapq", "percent_identity", "aln_block_len"],
                ascending=[True, False, False, False]
            ).iloc[0]

        results.append({
            "locus_id": locus["locus_id"],
            "contig": contig,
            "motifs_present": locus["motifs_present"],
            "strand_set": locus["strand_set"],
            "motif_count": locus["motif_count"],
            "contig_locus_start": locus_start,
            "contig_locus_stop": locus_stop,
            "block_contig_start": int(best["contig_start"]),
            "block_contig_stop": int(best["contig_stop"]),
            "ref_chr": best["ref_chr"],
            "chr_label": best["chr_label"],
            "block_ref_start": int(best["ref_start"]),
            "block_ref_stop": int(best["ref_stop"]),
            "block_strand": best["strand"],
            "block_mapq": int(best["mapq"]),
            "block_percent_identity": float(best["percent_identity"]),
            "block_aln_len": int(best["aln_block_len"]),
        })

    return pd.DataFrame(results)

def build_alignment_blocks(contigs_fasta: str | Path, ref_fasta: str | Path, out_dir: Path) -> pd.DataFrame:
    step("STEP 6A — Build alignment block table with minimap2 PAF")

    paf_file = out_dir / "contigs_vs_clib122.paf"
    run_minimap2_paf(contigs_fasta, ref_fasta, paf_file)

    raw_blocks = parse_paf_blocks(paf_file)
    raw_out = out_dir / "alignment_blocks_raw.csv"
    raw_blocks.to_csv(raw_out, index=False)
    info(f"Saved: {raw_out}")

    filt_blocks = filter_alignment_blocks(raw_blocks, min_block_len=5000, min_pid=90.0)
    filt_out = out_dir / "alignment_blocks_filtered.csv"
    filt_blocks.to_csv(filt_out, index=False)
    info(f"Saved: {filt_out}")

    if not filt_blocks.empty:
        print(filt_blocks.head(20).to_string(index=False))

    return filt_blocks

# =============================================================================
# STEP 6B - PAIR-SCORED, MULTI-FLANK BREAKPOINT PROJECTION
# =============================================================================

def parse_paf(paf_file: str | Path) -> pd.DataFrame:
    rows = []

    with Path(paf_file).open() as fh:
        for line in fh:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2]) + 1
            qend = int(parts[3])
            strand = parts[4]
            tname = parts[5]
            tlen = int(parts[6])
            tstart = int(parts[7]) + 1
            tend = int(parts[8])
            nmatch = int(parts[9])
            alnlen = int(parts[10])
            mapq = int(parts[11])

            pid = (nmatch / alnlen * 100) if alnlen > 0 else 0.0

            rows.append({
                "query": qname,
                "query_length": qlen,
                "query_start": qstart,
                "query_stop": qend,
                "strand": strand,
                "ref_chr": tname,
                "ref_length": tlen,
                "ref_start": tstart,
                "ref_stop": tend,
                "matches": nmatch,
                "aln_block_len": alnlen,
                "percent_identity": pid,
                "mapq": mapq,
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df["chr_label"] = df["ref_chr"].map(CHR_MAP).fillna(df["ref_chr"])
    return df


def extract_flank(seq: str, start: int, stop: int, flank_size: int = 10000) -> dict:
    """
    Contig coordinates are 1-based inclusive.
    Returns left and right flanks outside the cassette interval.
    """
    seq_len = len(seq)

    left_start = max(1, start - flank_size)
    left_stop = max(1, start - 1)

    right_start = min(seq_len, stop + 1)
    right_stop = min(seq_len, stop + flank_size)

    left_seq = seq[left_start - 1:left_stop] if left_stop >= left_start else ""
    right_seq = seq[right_start - 1:right_stop] if right_stop >= right_start else ""

    return {
        "left_start": left_start,
        "left_stop": left_stop,
        "left_seq": left_seq,
        "right_start": right_start,
        "right_stop": right_stop,
        "right_seq": right_seq,
    }


def filter_flank_candidates(
    paf_df: pd.DataFrame,
    expected_chr: str | None,
    expected_ref_start: int | None,
    expected_ref_stop: int | None,
    max_ref_distance_bp: int = 100000,
    min_pid: float = 90.0,
    min_aln_len: int = 500,
    top_n: int = 5,
) -> pd.DataFrame:
    """
    Keep top candidate alignments for one flank, constrained to the expected
    chromosome and local neighborhood of the assigned structural block.
    """
    if paf_df.empty:
        return pd.DataFrame()

    df = paf_df.copy()

    if expected_chr is not None:
        df = df[df["ref_chr"] == expected_chr].copy()

    if df.empty:
        return df

    df = df[
        (df["percent_identity"] >= min_pid) &
        (df["aln_block_len"] >= min_aln_len)
    ].copy()

    if df.empty:
        return df

    if expected_ref_start is not None and expected_ref_stop is not None:
        dists = []
        for _, r in df.iterrows():
            d = interval_distance(
                int(r["ref_start"]), int(r["ref_stop"]),
                int(expected_ref_start), int(expected_ref_stop)
            )
            dists.append(d)
        df["expected_block_distance_bp"] = dists

        near_df = df[df["expected_block_distance_bp"] <= max_ref_distance_bp].copy()
        if not near_df.empty:
            df = near_df

    df = df.sort_values(
        ["mapq", "percent_identity", "aln_block_len"],
        ascending=[False, False, False]
    ).head(top_n).reset_index(drop=True)

    return df


def score_flank_pair(
    left_row: pd.Series,
    right_row: pd.Series,
    expected_block_strand: str | None,
    max_anchor_separation_bp: int,
) -> Optional[dict]:
    """
    Score a left/right flank pair.
    Returns None if the pair is invalid.
    """
    if left_row["ref_chr"] != right_row["ref_chr"]:
        return None

    left_anchor = int(left_row["ref_stop"])
    right_anchor = int(right_row["ref_start"])

    interval_start = min(left_anchor, right_anchor)
    interval_stop = max(left_anchor, right_anchor)
    interval_width = interval_stop - interval_start + 1

    if interval_width > max_anchor_separation_bp:
        return None

    # Prefer strand consistency with structural block when available
    strand_pair = f"{left_row['strand']}/{right_row['strand']}"
    strand_penalty = 0.0
    if expected_block_strand is not None:
        if expected_block_strand == "+" and strand_pair != "+/+":
            strand_penalty = 10.0
        elif expected_block_strand == "-" and strand_pair != "-/-":
            strand_penalty = 10.0

    # Smaller interval is better. Higher identity and MAPQ are better.
    score = (
        float(left_row["mapq"]) +
        float(right_row["mapq"]) +
        float(left_row["percent_identity"]) +
        float(right_row["percent_identity"]) -
        (interval_width / 1000.0) -
        strand_penalty
    )

    return {
        "ref_chr": left_row["ref_chr"],
        "chr_label": left_row["chr_label"],
        "left_anchor": left_anchor,
        "right_anchor": right_anchor,
        "ref_start": interval_start,
        "ref_stop": interval_stop,
        "interval_width_bp": interval_width,
        "breakpoint_strand": strand_pair,
        "score": score,
        "left_mapq": int(left_row["mapq"]),
        "left_pid": float(left_row["percent_identity"]),
        "right_mapq": int(right_row["mapq"]),
        "right_pid": float(right_row["percent_identity"]),
    }


def choose_best_pair(
    left_candidates: pd.DataFrame,
    right_candidates: pd.DataFrame,
    expected_block_strand: str | None,
    max_anchor_separation_bp: int,
) -> Optional[dict]:
    """
    Evaluate all candidate left/right pairs and choose the best one.
    """
    if left_candidates.empty or right_candidates.empty:
        return None

    pair_results = []

    for _, lrow in left_candidates.iterrows():
        for _, rrow in right_candidates.iterrows():
            scored = score_flank_pair(
                lrow,
                rrow,
                expected_block_strand=expected_block_strand,
                max_anchor_separation_bp=max_anchor_separation_bp,
            )
            if scored is not None:
                pair_results.append(scored)

    if not pair_results:
        return None

    pair_df = pd.DataFrame(pair_results).sort_values(
        ["score", "interval_width_bp"],
        ascending=[False, True]
    )
    return pair_df.iloc[0].to_dict()


def choose_best_single_candidate(
    left_candidates: pd.DataFrame,
    right_candidates: pd.DataFrame,
) -> Optional[pd.Series]:
    """
    Fallback only when no valid pair exists.
    """
    candidates = []

    if not left_candidates.empty:
        candidates.append(left_candidates)
    if not right_candidates.empty:
        candidates.append(right_candidates)

    if not candidates:
        return None

    df = pd.concat(candidates, ignore_index=True)
    df = df.sort_values(
        ["mapq", "percent_identity", "aln_block_len"],
        ascending=[False, False, False]
    )
    return df.iloc[0]


def project_breakpoints_from_flanks(
    loci_df: pd.DataFrame,
    contigs: Dict[str, str],
    out_dir: Path,
    block_assignments_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Multi-pass breakpoint projection:
      - try flank sizes 10000, 5000, 2000
      - keep top candidate hits per flank
      - score left/right pairs
      - fallback to single flank only if needed
    """
    step("STEP 6B — Pair-scored, multi-flank breakpoint projection")

    if loci_df.empty:
        warn("No loci available for breakpoint projection.")
        return pd.DataFrame()

    results = []

    flank_sizes = [10000, 5000, 2000]
    max_anchor_separation_bp = 50000
    max_ref_distance_from_block_bp = 100000
    max_high_conf_interval_bp = 1000
    max_moderate_interval_bp = 5000

    for _, row in loci_df.iterrows():
        locus_id = row["locus_id"]
        contig = row["contig"]
        locus_start = int(row["locus_start"])
        locus_stop = int(row["locus_stop"])

        if contig not in contigs:
            warn(f"{locus_id}: contig not found in contig FASTA; skipping.")
            continue

        blk = block_assignments_df[block_assignments_df["locus_id"] == locus_id].copy()

        expected_chr = None
        expected_chr_label = None
        expected_ref_start = None
        expected_ref_stop = None
        expected_block_strand = None
        expected_block_pid = None
        expected_block_mapq = None

        if not blk.empty:
            blk_row = blk.iloc[0]
            expected_chr = blk_row["ref_chr"]
            expected_chr_label = blk_row["chr_label"]
            expected_ref_start = int(blk_row["block_ref_start"])
            expected_ref_stop = int(blk_row["block_ref_stop"])
            expected_block_strand = blk_row["block_strand"]
            expected_block_pid = float(blk_row["block_percent_identity"])
            expected_block_mapq = int(blk_row["block_mapq"])

        seq = contigs[contig]

        final_result = None
        final_status = None
        final_confidence = None
        final_flank_size = None
        left_best_mapq = None
        left_best_pid = None
        right_best_mapq = None
        right_best_pid = None

        for flank_size in flank_sizes:
            flank_info = extract_flank(seq, locus_start, locus_stop, flank_size=flank_size)

            left_candidates = pd.DataFrame()
            right_candidates = pd.DataFrame()

            if flank_info["left_seq"]:
                left_fasta = out_dir / f"{locus_id}_left_flank_{flank_size}.fasta"
                left_paf = out_dir / f"{locus_id}_left_flank_{flank_size}.paf"
                fasta_write({f"{locus_id}_LEFT_{flank_size}": flank_info["left_seq"]}, left_fasta)
                run_minimap2_paf(left_fasta, CONFIG["clib122_genome_fasta"], left_paf)
                left_df = parse_paf(left_paf)
                left_candidates = filter_flank_candidates(
                    left_df,
                    expected_chr=expected_chr,
                    expected_ref_start=expected_ref_start,
                    expected_ref_stop=expected_ref_stop,
                    max_ref_distance_bp=max_ref_distance_from_block_bp,
                    min_pid=90.0,
                    min_aln_len=500,
                    top_n=5,
                )

            if flank_info["right_seq"]:
                right_fasta = out_dir / f"{locus_id}_right_flank_{flank_size}.fasta"
                right_paf = out_dir / f"{locus_id}_right_flank_{flank_size}.paf"
                fasta_write({f"{locus_id}_RIGHT_{flank_size}": flank_info["right_seq"]}, right_fasta)
                run_minimap2_paf(right_fasta, CONFIG["clib122_genome_fasta"], right_paf)
                right_df = parse_paf(right_paf)
                right_candidates = filter_flank_candidates(
                    right_df,
                    expected_chr=expected_chr,
                    expected_ref_start=expected_ref_start,
                    expected_ref_stop=expected_ref_stop,
                    max_ref_distance_bp=max_ref_distance_from_block_bp,
                    min_pid=90.0,
                    min_aln_len=500,
                    top_n=5,
                )

            if not left_candidates.empty:
                left_best_mapq = int(left_candidates.iloc[0]["mapq"])
                left_best_pid = float(left_candidates.iloc[0]["percent_identity"])
            if not right_candidates.empty:
                right_best_mapq = int(right_candidates.iloc[0]["mapq"])
                right_best_pid = float(right_candidates.iloc[0]["percent_identity"])

            best_pair = choose_best_pair(
                left_candidates,
                right_candidates,
                expected_block_strand=expected_block_strand,
                max_anchor_separation_bp=max_anchor_separation_bp,
            )

            if best_pair is not None:
                final_result = best_pair
                final_flank_size = flank_size

                if best_pair["interval_width_bp"] <= max_high_conf_interval_bp:
                    final_status = "paired_flanks"
                    final_confidence = "high"
                elif best_pair["interval_width_bp"] <= max_moderate_interval_bp:
                    final_status = "paired_flanks_broad_interval"
                    final_confidence = "moderate"
                else:
                    final_status = "paired_flanks_broad_interval"
                    final_confidence = "moderate"

                # stop early on strongest solution
                if final_confidence == "high":
                    break

        # fallback if no valid pair found after all flank sizes
        if final_result is None:
            best_single = choose_best_single_candidate(left_candidates, right_candidates)

            if best_single is None:
                warn(f"{locus_id}: no usable flank alignments found.")
                results.append({
                    "locus_id": locus_id,
                    "contig": contig,
                    "motifs_present": row["motifs_present"],
                    "strand_set": row["strand_set"],
                    "motif_count": row["motif_count"],
                    "contig_locus_start": locus_start,
                    "contig_locus_stop": locus_stop,
                    "expected_block_ref_chr": expected_chr,
                    "expected_block_chr_label": expected_chr_label,
                    "expected_block_ref_start": expected_ref_start,
                    "expected_block_ref_stop": expected_ref_stop,
                    "expected_block_strand": expected_block_strand,
                    "expected_block_mapq": expected_block_mapq,
                    "expected_block_pid": expected_block_pid,
                    "ref_chr": expected_chr,
                    "chr_label": expected_chr_label,
                    "ref_start": pd.NA,
                    "ref_stop": pd.NA,
                    "breakpoint_strand": pd.NA,
                    "projection_status": "unresolved",
                    "projection_confidence": "low",
                    "projection_interval_width_bp": pd.NA,
                    "selected_flank_size_bp": pd.NA,
                    "left_flank_mapq": left_best_mapq,
                    "left_flank_pid": left_best_pid,
                    "right_flank_mapq": right_best_mapq,
                    "right_flank_pid": right_best_pid,
                })
                continue

            final_result = {
                "ref_chr": best_single["ref_chr"],
                "chr_label": best_single["chr_label"],
                "ref_start": int(best_single["ref_start"]),
                "ref_stop": int(best_single["ref_stop"]),
                "breakpoint_strand": str(best_single["strand"]),
                "interval_width_bp": int(best_single["ref_stop"]) - int(best_single["ref_start"]) + 1,
            }
            final_status = "single_flank_fallback"
            final_confidence = "low"
            final_flank_size = flank_sizes[-1]

        results.append({
            "locus_id": locus_id,
            "contig": contig,
            "motifs_present": row["motifs_present"],
            "strand_set": row["strand_set"],
            "motif_count": row["motif_count"],

            "contig_locus_start": locus_start,
            "contig_locus_stop": locus_stop,

            "expected_block_ref_chr": expected_chr,
            "expected_block_chr_label": expected_chr_label,
            "expected_block_ref_start": expected_ref_start,
            "expected_block_ref_stop": expected_ref_stop,
            "expected_block_strand": expected_block_strand,
            "expected_block_mapq": expected_block_mapq,
            "expected_block_pid": expected_block_pid,

            "ref_chr": final_result["ref_chr"],
            "chr_label": final_result["chr_label"],
            "ref_start": final_result["ref_start"],
            "ref_stop": final_result["ref_stop"],
            "breakpoint_strand": final_result["breakpoint_strand"],
            "projection_status": final_status,
            "projection_confidence": final_confidence,
            "projection_interval_width_bp": final_result["interval_width_bp"],
            "selected_flank_size_bp": final_flank_size,

            "left_flank_mapq": left_best_mapq,
            "left_flank_pid": left_best_pid,
            "right_flank_mapq": right_best_mapq,
            "right_flank_pid": right_best_pid,
        })

    proj_df = pd.DataFrame(results)
    out = out_dir / "projected_loci_breakpoint_based.csv"
    proj_df.to_csv(out, index=False)
    info(f"Saved: {out}")

    if not proj_df.empty:
        print(proj_df.to_string(index=False))

    return proj_df

# =============================================================================
# STEP 7 - PARSE CLIB122 GFF3
# =============================================================================

def parse_gff3_to_table(gff_path: str | Path, out_dir: Path) -> pd.DataFrame:
    step("STEP 7 — Parse CLIB122 GFF3 into a clean annotation table")

    rows = []
    with Path(gff_path).open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attrs = parts
            attr_dict = {}
            for field in attrs.split(";"):
                if "=" in field:
                    k, v = field.split("=", 1)
                    attr_dict[k] = v

            rows.append({
                "chr": seqid,
                "source": source,
                "Type": feature_type,
                "start": int(start),
                "stop": int(end),
                "Direction": strand,
                "Name": attr_dict.get("Name", ""),
                "gene": attr_dict.get("gene", ""),
                "product": attr_dict.get("product", ""),
                "note": attr_dict.get("Note", ""),
                "dbxref": attr_dict.get("Dbxref", ""),
                "locus_tag": attr_dict.get("locus_tag") or attr_dict.get("ID") or "",
                "attributes": attrs,
            })

    ann = pd.DataFrame(rows)
    if not ann.empty:
        ann["chr_label"] = ann["chr"].map(CHR_MAP).fillna(ann["chr"])
    out = out_dir / "CLIB122_annotation_table.csv"
    ann.to_csv(out, index=False)
    info(f"Saved: {out}")
    info(f"Parsed {len(ann)} annotation rows")
    return ann


# =============================================================================
# STEP 8 - BUILD MAPPED ESSENTIAL TABLE
# =============================================================================

def build_mapped_essential_table(out_dir: Path) -> pd.DataFrame:
    step("STEP 8 — Build mapped essential-gene table (CLIB89 -> CLIB122)")

    essential = pd.read_excel(CONFIG["clib89_essential_table"], skiprows=1)
    mapping = read_table(CONFIG["clib89_clib122_mapping_table"])

    essential.columns = essential.columns.astype(str).str.strip()
    mapping.columns = mapping.columns.astype(str).str.strip()

    essential = essential[essential["Gen80Classification"] == "Essential"].copy()
    essential_id_col = "Name"
    if essential_id_col not in essential.columns:
        fail(f"Expected column '{essential_id_col}' not found in essential table.")

    yali1_map_col = detect_headerless_ids(mapping, ["YALI1 Locus ID"])
    yali0_map_col = detect_headerless_ids(mapping, ["YALI0 locus ID"])
    if yali1_map_col is None or yali0_map_col is None:
        for c in mapping.columns:
            if "yali1" in c.lower():
                yali1_map_col = c
            if "yali0" in c.lower():
                yali0_map_col = c
    if yali1_map_col is None or yali0_map_col is None:
        fail("Could not identify the YALI1 and YALI0 columns in the mapping table.")

    essential["YALI1_norm"] = essential[essential_id_col].astype(str).str.strip()
    mapping["YALI1_norm"] = mapping[yali1_map_col].astype(str).str.strip()
    mapping["YALI0_norm"] = mapping[yali0_map_col].astype(str).str.strip().apply(normalize_yali0)

    merged = essential.merge(mapping, on="YALI1_norm", how="left", suffixes=("_essential", "_mapping"))
    out_full = out_dir / "mapped_essential_genes_full.csv"
    merged.to_csv(out_full, index=False)

    essential_only = merged[
        merged["YALI0_norm"].notna() &
        (merged["YALI0_norm"].astype(str).str.strip() != "") &
        (merged["YALI0_norm"].astype(str).str.strip().str.lower() != "nan")
    ].copy()

    out_simple = out_dir / "mapped_essential_gene_list.csv"
    essential_only[["YALI1_norm", "YALI0_norm"]].drop_duplicates().to_csv(out_simple, index=False)

    info(f"Saved: {out_full}")
    info(f"Saved: {out_simple}")
    info(f"Mapped essential rows: {len(essential_only)}")
    return essential_only


# =============================================================================
# STEP 9 - EXTRACT INSTABILITY FEATURES
# =============================================================================

def extract_instability_features(annotation_df: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    step("STEP 9 — Extract instability-associated features from CLIB122 annotation")

    ann = annotation_df.copy()
    instability = ann[
        (ann["Type"].isin(["mobile_element", "repeat_region", "tRNA", "rRNA"])) |
        (ann["Name"].astype(str).str.contains("Centromere|repeat|LTR", case=False, na=False)) |
        (ann["product"].astype(str).str.contains("repeat|LTR|retrotransposon", case=False, na=False)) |
        (ann["note"].astype(str).str.contains("repeat|LTR|retrotransposon", case=False, na=False))
    ].copy()

    def classify(row):
        t = str(row["Type"])
        n = str(row.get("Name", ""))
        p = str(row.get("product", ""))
        note = str(row.get("note", ""))
        text = " ".join([t, n, p, note]).lower()

        if "mobile_element" in t:
            return "high_mobile"
        if "repeat_region" in t or "repeat" in text or "ltr" in text or "retrotransposon" in text:
            return "high_repeat"
        if "centromere" in text:
            return "structural_centromere"
        if t == "tRNA":
            return "moderate_tRNA"
        if t == "rRNA":
            return "moderate_rRNA"
        return "other"

    instability["instability_class"] = instability.apply(classify, axis=1)
    if not instability.empty:
        instability["chr_label"] = instability["chr"].map(CHR_MAP).fillna(instability["chr"])

    keep_cols = [c for c in [
        "chr", "chr_label", "source", "Type", "start", "stop", "Direction",
        "Name", "gene", "product", "note", "dbxref", "locus_tag",
        "instability_class", "attributes"
    ] if c in instability.columns]
    instability = instability[keep_cols].copy()

    out = out_dir / "instability_features.csv"
    instability.to_csv(out, index=False)
    info(f"Saved: {out}")
    info(f"Instability features: {len(instability)}")
    return instability


# =============================================================================
# WINDOW HELPER
# =============================================================================

def features_in_window(query_chr, query_start, query_stop, annotation_df, window_bp):
    sub = annotation_df[annotation_df["chr"] == query_chr].copy()
    if sub.empty:
        return sub

    sub = sub[
        (sub["stop"] >= query_start - window_bp) &
        (sub["start"] <= query_stop + window_bp)
    ].copy()
    if sub.empty:
        return sub

    sub["distance_bp"] = [
        interval_distance(query_start, query_stop, int(s), int(e))
        for s, e in zip(sub["start"], sub["stop"])
    ]
    sub["overlaps_locus"] = sub["distance_bp"] == 0
    return sub.sort_values(["distance_bp", "start"]).copy()


# =============================================================================
# STEP 10 - CONFIDENCE-AWARE NEIGHBORHOOD ANALYSIS
# =============================================================================

def format_feature_list(df, label_col, distance_col="distance_bp", max_items=6):
    if df is None or df.empty:
        return ""

    items = []
    for _, row in df.head(max_items).iterrows():
        label = str(row.get(label_col, "")).strip()
        if not label:
            label = str(row.get("product", "")).strip()
        if not label:
            label = str(row.get("Name", "")).strip()
        if not label:
            label = "feature"

        dist = row.get(distance_col, pd.NA)
        if pd.notna(dist):
            items.append(f"{label} ({int(dist)} bp)")
        else:
            items.append(label)

    return "; ".join(items)


def analyze_locus_neighborhoods(
    projected_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    mapped_essential_df: pd.DataFrame,
    instability_df: pd.DataFrame,
    out_dir: Path,
) -> None:
    step("STEP 10 — Confidence-aware neighborhood analysis")

    if projected_df.empty:
        warn("No projected loci available for neighborhood analysis.")
        return

    gene_window_bp = 20000
    instability_window_bp = 50000

    ann2 = annotation_df.copy()
    ann2["locus_tag_norm"] = ann2["locus_tag"].astype(str).str.strip().apply(normalize_yali0)

    # Collapse CDS rows to one gene-level row per locus_tag
    cds_like = ann2[ann2["Type"] == "CDS"].copy()

    def first_nonempty(series):
        for x in series:
            if pd.notna(x) and str(x).strip() != "":
                return x
        return ""

    cds_gene_level = (
        cds_like
        .groupby("locus_tag", dropna=False)
        .agg({
            "chr": "first",
            "source": "first",
            "Direction": "first",
            "Name": first_nonempty,
            "gene": first_nonempty,
            "product": first_nonempty,
            "note": first_nonempty,
            "dbxref": first_nonempty,
            "start": "min",
            "stop": "max",
        })
        .reset_index()
    )

    cds_gene_level["Type"] = "CDS_gene_collapsed"
    cds_gene_level["locus_tag_norm"] = cds_gene_level["locus_tag"].astype(str).str.strip().apply(normalize_yali0)

    essential_set = set(
        mapped_essential_df["YALI0_norm"].astype(str).str.strip().apply(normalize_yali0)
    )
    cds_gene_level["is_essential"] = cds_gene_level["locus_tag_norm"].isin(essential_set)

    cds_window_rows = []
    instability_window_rows = []
    summary_rows = []

    for _, row in projected_df.iterrows():
        locus_id = row["locus_id"]
        ref_chr = row["ref_chr"]
        chr_label = row.get("chr_label", CHR_MAP.get(ref_chr, ref_chr))

        if pd.isna(ref_chr) or pd.isna(row.get("ref_start", pd.NA)) or pd.isna(row.get("ref_stop", pd.NA)):
            continue

        qstart = int(row["ref_start"])
        qstop = int(row["ref_stop"])

        projection_status = str(row.get("projection_status", ""))
        projection_confidence = str(row.get("projection_confidence", ""))
        projection_interval_width_bp = int(row.get("projection_interval_width_bp", qstop - qstart + 1))

        allow_exact_overlap = projection_confidence == "high"

        cds_window = features_in_window(
            ref_chr, qstart, qstop, cds_gene_level, gene_window_bp
        )

        if not cds_window.empty:
            # Only high-confidence projections are allowed to claim exact overlap
            if not allow_exact_overlap:
                cds_window["overlaps_locus"] = False

            interval_hits = cds_window[
                (cds_window["stop"] >= qstart) & (cds_window["start"] <= qstop)
            ].copy()

            flanking = cds_window[cds_window["distance_bp"] > 0].copy()
            essential = cds_window[cds_window["is_essential"]].copy()

            for _, g in cds_window.iterrows():
                cds_window_rows.append({
                    "locus_id": locus_id,
                    "contig": row["contig"],
                    "ref_chr": ref_chr,
                    "chr_label": chr_label,
                    "locus_start": qstart,
                    "locus_stop": qstop,
                    "motifs_present": row.get("motifs_present", ""),
                    "projection_status": projection_status,
                    "projection_confidence": projection_confidence,
                    "projection_interval_width_bp": projection_interval_width_bp,
                    "locus_tag": g["locus_tag"],
                    "product": g.get("product", ""),
                    "note": g.get("note", ""),
                    "cds_start": int(g["start"]),
                    "cds_stop": int(g["stop"]),
                    "distance_bp": int(g["distance_bp"]),
                    "overlaps_locus": bool(g["overlaps_locus"]),
                    "within_projected_interval": bool((int(g["stop"]) >= qstart) and (int(g["start"]) <= qstop)),
                    "is_essential": bool(g["is_essential"]),
                })
        else:
            interval_hits = pd.DataFrame()
            flanking = pd.DataFrame()
            essential = pd.DataFrame()

        inst_window = features_in_window(
            ref_chr, qstart, qstop, instability_df, instability_window_bp
        )

        inst_close = inst_window[inst_window["distance_bp"] <= 5000].copy() if not inst_window.empty else pd.DataFrame()
        inst_local = inst_window[inst_window["distance_bp"] <= 20000].copy() if not inst_window.empty else pd.DataFrame()

        if not inst_window.empty:
            for _, g in inst_window.iterrows():
                instability_window_rows.append({
                    "locus_id": locus_id,
                    "contig": row["contig"],
                    "ref_chr": ref_chr,
                    "chr_label": chr_label,
                    "locus_start": qstart,
                    "locus_stop": qstop,
                    "projection_status": projection_status,
                    "projection_confidence": projection_confidence,
                    "projection_interval_width_bp": projection_interval_width_bp,
                    "feature_type": g["Type"],
                    "product": g.get("product", ""),
                    "note": g.get("note", ""),
                    "instability_class": g["instability_class"],
                    "feature_start": int(g["start"]),
                    "feature_stop": int(g["stop"]),
                    "distance_bp": int(g["distance_bp"]),
                    "overlaps_locus": bool(g["overlaps_locus"]) if allow_exact_overlap else False,
                    "within_projected_interval": bool((int(g["stop"]) >= qstart) and (int(g["start"]) <= qstop)),
                })

        exact_overlap_str = ""
        if allow_exact_overlap:
            exact_overlap_df = cds_window[cds_window["overlaps_locus"]] if not cds_window.empty else pd.DataFrame()
            exact_overlap_str = format_feature_list(exact_overlap_df, "locus_tag")

        genes_in_interval_str = format_feature_list(interval_hits, "locus_tag")
        flanking_str = format_feature_list(flanking, "locus_tag")
        essential_str = format_feature_list(essential, "locus_tag")
        inst_close_str = format_feature_list(inst_close, "product") if not inst_close.empty else ""
        inst_local_str = format_feature_list(inst_local, "product") if not inst_local.empty else ""
        num_genes_in_interval = len(interval_hits) if interval_hits is not None else 0
        num_flanking_genes = len(flanking) if flanking is not None else 0

        if projection_confidence == "high":
            interpretation = "Precise insertion site (high confidence)"
        elif projection_confidence == "moderate":
            interpretation = "Localized insertion region (~5–10 kb)"
        else:
            interpretation = "Chromosomal assignment only (low confidence)"

        summary_rows.append({
            "locus_id": locus_id,
            "contig": row["contig"],
            "chr_label": chr_label,
            "locus_start": qstart,
            "locus_stop": qstop,
            "motifs_present": row.get("motifs_present", ""),

            "projection_status": projection_status,
            "projection_confidence": projection_confidence,
            "projection_interval_width_bp": projection_interval_width_bp,

            "interpretation_level": interpretation,

            "exact_overlapping_cds_high_confidence_only": exact_overlap_str,
            "genes_within_insertion_region": genes_in_interval_str,
            "nearest_flanking_genes": flanking_str,
            "essential_genes_in_window": essential_str,

            "instability_features_immediate": inst_close_str,
            "instability_features_local": inst_local_str,

            "num_genes_in_interval": num_genes_in_interval,
            "num_flanking_genes": num_flanking_genes,
            "num_instability_features_immediate": len(inst_close),
            "num_instability_features_local": len(inst_local),
        })

    cds_window_df = pd.DataFrame(cds_window_rows)
    instability_window_df = pd.DataFrame(instability_window_rows)
    summary_df = pd.DataFrame(summary_rows)

    cds_window_df.to_csv(out_dir / "locus_cds_window.csv", index=False)
    instability_window_df.to_csv(out_dir / "locus_instability_window.csv", index=False)
    summary_df.to_csv(out_dir / "locus_summary_compact.csv", index=False)

    info(f"Saved: {out_dir / 'locus_cds_window.csv'}")
    info(f"Saved: {out_dir / 'locus_instability_window.csv'}")
    info(f"Saved: {out_dir / 'locus_summary_compact.csv'}")

# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("\nYarrowia engineered-strain genomic context pipeline")
    print("Single-file, step-by-step workflow\n")

    out_dir = ensure_dir(CONFIG["output_dir"])
    pd.Series({k: str(v) for k, v in CONFIG.items()}).to_csv(out_dir / "run_config.csv", header=["value"])

    motif_fasta, motifs = build_motif_fasta(out_dir)
    contigs, refs = load_inputs()
    assignment_df = assign_contigs_to_reference(contigs, refs, out_dir)
    hits_df = detect_motifs_on_contigs(contigs, motifs, out_dir)
    loci_df = group_hits_into_loci(hits_df, out_dir)

    blocks_df = build_alignment_blocks(
        CONFIG["contigs_fasta"],
        CONFIG["clib122_genome_fasta"],
        out_dir
    )

    locus_blocks_df = assign_loci_to_blocks(loci_df, blocks_df)
    locus_blocks_df.to_csv(out_dir / "loci_assigned_blocks.csv", index=False)
    info(f"Saved: {out_dir / 'loci_assigned_blocks.csv'}")


    projected_df = project_breakpoints_from_flanks(
    loci_df,
    contigs,
    out_dir,
    block_assignments_df=locus_blocks_df
    )

    ann_df = parse_gff3_to_table(CONFIG["clib122_gff3"], out_dir)
    mapped_essential_df = build_mapped_essential_table(out_dir)
    instability_df = extract_instability_features(ann_df, out_dir)
    analyze_locus_neighborhoods(projected_df, ann_df, mapped_essential_df, instability_df, out_dir)

    step("STEP 11 — Final interpretation note")
    print(
        "Alignment blocks provide structural context (chromosome, orientation, inversions),\n"
        "while breakpoint-based projections provide a conservative inferred insertion interval.\n"
        "For broad or low-confidence projected intervals, treat overlap calls as unresolved and\n"
        "interpret the outputs as local neighborhood candidates rather than definitive disruptions.\n"
    )

    step("DONE")
    print("Main output files:")
    for fname in [
        "contig_assignment.csv",
        "motif_hits.csv",
        "integration_loci.csv",
        "alignment_blocks_raw.csv",
        "alignment_blocks_filtered.csv",
        "loci_assigned_blocks.csv",
        "projected_loci_breakpoint_based.csv",
        "CLIB122_annotation_table.csv",
        "mapped_essential_genes_full.csv",
        "mapped_essential_gene_list.csv",
        "instability_features.csv",
        "locus_cds_window.csv",
        "locus_instability_window.csv",
        "locus_summary_compact.csv",
    ]:
        print(" -", out_dir / fname)

if __name__ == "__main__":
    main()
