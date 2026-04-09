#!/usr/bin/env python3
"""
Generate synthetic test BAM files and BED files for cen-tss-analysis pipeline testing.

Usage:
    python tests/make_test_data.py

Outputs written to tests/data/:
    centromere.bed          Centromere regions BED  (step 1 input)
    transposon.bed          Transposon regions BED  (step 5 input)
    step1_input.bam         STAR-like genome BAM    (step 1 input)
    step3_input.bam         Bowtie2-like remap BAM  (step 3 / pyselectal.py input)
    step4_input.bam         TSS-selected BAM        (step 4 / xm_filter input)
    step5_input.bam         XM-filtered BAM         (step 5 / transposon filter input)

Each BAM is coordinate-sorted and indexed (.bai).

Design notes
------------
step1_input.bam  uses chr1 (genome reference, length 100 000).
steps 3-5 BAMs   use HOR_chr1 (centromeric reference, length 100 000).

Expected outcomes after each step (used by run_tests.sh):
    step 1  → 8 reads in output FASTQ  (5 centromeric + 3 unmapped)
    step 3  → 3 reads in output BAM    (2 fwd + 1 rev with exactly 1 5'-soft-clipped G)
    step 4  → 2 reads in output BAM    (XM/readlen <= 0.10)
    step 5  → 3 reads in output BAM    (not overlapping transposon HOR_chr1:1000-2000)

Coordinate conventions
----------------------
All BAM positions are 0-based (pysam reference_start).
BED files use 0-based half-open intervals [start, end).
bedtools intersect -abam: uses mapped region only (soft-clips excluded).

For a read with CIGAR 1S+49M at reference_start=X:
    mapped region = [X, X+49)
    soft-clip base = X-1 on the reference (not counted by bedtools)
"""

import os
import sys
import pysam

OUT_DIR = os.path.join(os.path.dirname(__file__), "data")

# ── CIGAR operation codes (pysam) ─────────────────────────────────────────────
MATCH = 0  # M
SOFT  = 4  # S


# ── helpers ───────────────────────────────────────────────────────────────────

def ensure_out_dir():
    os.makedirs(OUT_DIR, exist_ok=True)


def write_bam(filename, header_dict, build_fn):
    """
    Write a coordinate-sorted, indexed BAM.

    build_fn(bam) is called with an open AlignmentFile in 'wb' mode; it should
    call bam.write() for every record.  Sorting is done via pysam.sort on a
    temporary unsorted file.
    """
    path = os.path.join(OUT_DIR, filename)
    tmp  = path + ".unsorted.bam"

    with pysam.AlignmentFile(tmp, "wb", header=header_dict) as bam:
        build_fn(bam)

    pysam.sort("-o", path, tmp)
    os.remove(tmp)
    pysam.index(path)
    print(f"  wrote {path}")
    return path


def aln(bam, name, ref_id, pos, cigar, seq, *, is_reverse=False,
        is_unmapped=False, mapq=60, tags=None):
    """
    Build and return one pysam.AlignedSegment attached to *bam*.

    Parameters
    ----------
    bam         : open pysam.AlignmentFile
    name        : query name (str)
    ref_id      : 0-based reference index (-1 for unmapped)
    pos         : 0-based leftmost mapped position
    cigar       : list of (op, length) tuples, or None for unmapped
    seq         : query sequence string (full, including soft-clipped bases)
    is_reverse  : set FLAG 0x10
    is_unmapped : set FLAG 0x4; ignores ref_id/pos/cigar
    mapq        : mapping quality
    tags        : list of (tag, value) tuples, e.g. [("XM", 2), ("AS", 40)]
    """
    a = pysam.AlignedSegment(bam.header)
    a.query_name      = name
    a.query_sequence  = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))

    if is_unmapped:
        a.flag             = 4
        a.reference_id     = -1
        a.reference_start  = 0
        a.mapping_quality  = 0
        a.cigar            = None
    else:
        a.flag            = 16 if is_reverse else 0
        a.reference_id    = ref_id
        a.reference_start = pos
        a.mapping_quality = mapq
        a.cigar           = cigar

    if tags:
        a.set_tags(tags)

    return a


# ── BED files ─────────────────────────────────────────────────────────────────

def make_centromere_bed():
    path = os.path.join(OUT_DIR, "centromere.bed")
    with open(path, "w") as f:
        # One centromere band on chr1, 500–2000
        f.write("chr1\t500\t2000\t.\t.\t+\n")
    print(f"  wrote {path}")


def make_transposon_bed():
    path = os.path.join(OUT_DIR, "transposon.bed")
    with open(path, "w") as f:
        # One transposon on the centromeric reference, 1000–2000
        f.write("HOR_chr1\t1000\t2000\n")
    print(f"  wrote {path}")


# ── Step 1 input BAM ──────────────────────────────────────────────────────────
# Simulates STAR-mapped full-genome BAM (reference: chr1).
#
# Contents:
#   5 mapped reads at chr1:900-1090  → inside centromere.bed (chr1:500-2000)
#   3 unmapped reads                 → captured via samtools -f 4
#   4 mapped reads at chr1:3000-3150 → outside centromere, should be excluded
#
# Expected step-1 output: FASTQ with 5 + 3 = 8 reads.

def make_step1_input():
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 100000}],
        "PG": [{"ID": "STAR", "PN": "STAR", "VN": "2.7.10a"}],
    }

    def build(bam):
        # 5 mapped reads inside centromere (chr1:900-…)
        for i in range(5):
            bam.write(aln(bam,
                name=f"cen_read_{i+1}",
                ref_id=0, pos=900 + i * 38,
                cigar=[(MATCH, 50)],
                seq="A" * 50,
                tags=[("NH", 1), ("HI", 1), ("nM", 0)]))

        # 3 unmapped reads
        for i in range(3):
            bam.write(aln(bam,
                name=f"unmapped_read_{i+1}",
                ref_id=-1, pos=0,
                cigar=None,
                seq="T" * 40,
                is_unmapped=True))

        # 4 mapped reads outside centromere (chr1:3000-…)
        for i in range(4):
            bam.write(aln(bam,
                name=f"non_cen_read_{i+1}",
                ref_id=0, pos=3000 + i * 50,
                cigar=[(MATCH, 50)],
                seq="C" * 50,
                tags=[("NH", 1), ("HI", 1), ("nM", 0)]))

    write_bam("step1_input.bam", header, build)


# ── Step 3 input BAM ──────────────────────────────────────────────────────────
# Simulates bowtie2-remapped BAM to HOR_chr1.
# Tests pyselectal.py -n 1 -m 1 -x G  (exactly 1 5'-soft-clipped G).
#
# PASS criteria:
#   Forward: CIGAR starts with (SOFT,1), first base of query seq == 'G'
#   Reverse: CIGAR ends   with (SOFT,1), last  base of query seq == 'C'
#            (the soft-clipped base on the reverse strand is the reverse
#            complement of G, i.e. C; pysam stores the read sequence in its
#            original orientation before reverse-complementing)
#
# read name          expected  reason
# -----------------  --------  ------------------------------------------------
# fwd_pass_1G        PASS      1S+50M, first base G  — canonical CAGE fwd read
# fwd_pass_2G        PASS      1S+50M, first base G  — second fwd example
# rev_pass_1C        PASS      50M+1S, last base C   — canonical CAGE rev read
# fwd_fail_notG      FAIL      1S+50M but first base is A, not G
# fwd_fail_2Sclip    FAIL      2S+49M — two-base soft-clip; -n 1 -m 1 requires exactly 1
# fwd_fail_noSclip   FAIL      51M — no soft-clip at all
# rev_fail_notC      FAIL      50M+1S but last base is A, not C

def make_step3_input():
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "HOR_chr1", "LN": 100000}],
        "PG": [{"ID": "bowtie2", "PN": "bowtie2", "VN": "2.5.1"}],
    }

    def build(bam):
        # ── PASS ──────────────────────────────────────────────────────────────
        bam.write(aln(bam, "fwd_pass_1G", 0, 500,
            cigar=[(SOFT, 1), (MATCH, 50)],
            seq="G" + "A" * 50,
            tags=[("XM", 0), ("AS", 45)]))

        bam.write(aln(bam, "fwd_pass_2G", 0, 600,
            cigar=[(SOFT, 1), (MATCH, 50)],
            seq="G" + "T" * 50,
            tags=[("XM", 0), ("AS", 44)]))

        # Reverse: CIGAR is 50M+1S; last base of stored seq = C
        bam.write(aln(bam, "rev_pass_1C", 0, 700,
            cigar=[(MATCH, 50), (SOFT, 1)],
            seq="A" * 50 + "C",
            is_reverse=True,
            tags=[("XM", 0), ("AS", 43)]))

        # ── FAIL ──────────────────────────────────────────────────────────────
        # 1S but first base is A (not G)
        bam.write(aln(bam, "fwd_fail_notG", 0, 800,
            cigar=[(SOFT, 1), (MATCH, 50)],
            seq="A" + "T" * 50,
            tags=[("XM", 0), ("AS", 40)]))

        # 2S + 49M (wrong soft-clip length)
        bam.write(aln(bam, "fwd_fail_2Sclip", 0, 900,
            cigar=[(SOFT, 2), (MATCH, 49)],
            seq="GG" + "A" * 49,
            tags=[("XM", 0), ("AS", 39)]))

        # No soft-clip (51M)
        bam.write(aln(bam, "fwd_fail_noSclip", 0, 1000,
            cigar=[(MATCH, 51)],
            seq="G" + "A" * 50,
            tags=[("XM", 0), ("AS", 48)]))

        # Reverse: last base is A, not C
        bam.write(aln(bam, "rev_fail_notC", 0, 1100,
            cigar=[(MATCH, 50), (SOFT, 1)],
            seq="A" * 50 + "A",
            is_reverse=True,
            tags=[("XM", 0), ("AS", 38)]))

    write_bam("step3_input.bam", header, build)


# ── Step 4 input BAM ──────────────────────────────────────────────────────────
# Tests xm_filter.sh: keep reads where XM / read_length <= 0.10.
# All reads are 30 bp with CIGAR 1S+29M (already TSS-selected).
#
# xm_filter.sh computes: mismatch_fraction = XM / len(SEQ)
# Threshold: fraction must be <= (1 - min_match_frac) = 1 - 0.90 = 0.10
#
# read name             expected  reason
# --------------------  --------  ----------------------------------------------
# xm_pass_0mm           PASS      XM=0  (0/30  = 0.000, no mismatches)
# xm_pass_3mm_boundary  PASS      XM=3  (3/30  = 0.100, exactly at 10% limit)
# xm_fail_4mm           FAIL      XM=4  (4/30  = 0.133, just over 10%)
# xm_fail_10mm          FAIL      XM=10 (10/30 = 0.333, well over 10%)
# xm_drop_noXM          FAIL      no XM tag — awk never matches; read dropped
#
# Note: xm_drop_noXM is silently dropped (not explicitly failed). Reads without
# XM tag arise when a mapper does not output XM (e.g. some bowtie2 modes). The
# pipeline currently treats them as non-passing, same as high-mismatch reads.

def make_step4_input():
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "HOR_chr1", "LN": 100000}],
        "PG": [{"ID": "bowtie2", "PN": "bowtie2", "VN": "2.5.1"}],
    }
    SEQ30 = "G" + "A" * 29   # 30 bp; starts with G (already TSS-selected)

    def build(bam):
        bam.write(aln(bam, "xm_pass_0mm", 0, 500,
            cigar=[(SOFT, 1), (MATCH, 29)], seq=SEQ30,
            tags=[("XM", 0), ("AS", 55)]))

        bam.write(aln(bam, "xm_pass_3mm_boundary", 0, 550,
            cigar=[(SOFT, 1), (MATCH, 29)], seq=SEQ30,
            tags=[("XM", 3), ("AS", 50)]))

        bam.write(aln(bam, "xm_fail_4mm", 0, 600,
            cigar=[(SOFT, 1), (MATCH, 29)], seq=SEQ30,
            tags=[("XM", 4), ("AS", 46)]))

        bam.write(aln(bam, "xm_fail_10mm", 0, 650,
            cigar=[(SOFT, 1), (MATCH, 29)], seq=SEQ30,
            tags=[("XM", 10), ("AS", 35)]))

        # No XM tag — awk loop never matches /^XM:i:/, read is silently dropped
        bam.write(aln(bam, "xm_drop_noXM", 0, 700,
            cigar=[(SOFT, 1), (MATCH, 29)], seq=SEQ30,
            tags=[("AS", 55)]))

    write_bam("step4_input.bam", header, build)


# ── Step 5 input BAM ──────────────────────────────────────────────────────────
# Tests filter_transposons.sh: remove reads overlapping transposon.bed
# (HOR_chr1:1000-2000).  All reads are 50 bp with CIGAR 1S+49M.
#
# Coordinate arithmetic for CIGAR 1S+49M at reference_start=X:
#   soft-clip occupies one query base but NO reference position
#   mapped region = [X, X+49)  (0-based half-open, as bedtools sees it)
#
# transposon BED: HOR_chr1  1000  2000  → interval [1000, 2000) 0-based
#
# read name               expected  strand  mapped interval  overlap with [1000,2000)?
# ----------------------  --------  ------  ---------------  -------------------------
# transp_pass_upstream    PASS      fwd     [500,  549)      no  — entirely before 1000
# transp_pass_downstream  PASS      fwd     [2500, 2549)     no  — entirely after 2000
# transp_pass_boundary    PASS      rev     [2000, 2049)     no  — starts exactly at 2000
# transp_fail_inside      FAIL      fwd     [1100, 1149)     yes — fully inside
# transp_fail_overlap_start FAIL    fwd     [970,  1019)     yes — overlaps start [970,1000)∩[1000,2000)
# transp_fail_overlap_end FAIL      fwd     [1970, 2019)     yes — overlaps end   [1970,2000)
#
# transp_pass_boundary uses is_reverse=True so that step 5 output contains both
# forward and reverse reads — required for step 6 to produce non-empty .fwd.bw
# and .rev.bw bigWigs.

def make_step5_input():
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "HOR_chr1", "LN": 100000}],
        "PG": [{"ID": "bowtie2", "PN": "bowtie2", "VN": "2.5.1"}],
    }
    # All reads: 50 bp, CIGAR 1S+49M (already TSS-selected).
    # Forward: seq = G + A*49 (5'-G soft-clip)
    # Reverse: seq = A*49 + C (5'-G soft-clip stored as C on reverse strand)
    SEQ_FWD = "G" + "A" * 49
    SEQ_REV = "A" * 49 + "C"

    def build(bam):
        # PASS: mapped [500, 549) — entirely before transposon [1000, 2000)
        bam.write(aln(bam, "transp_pass_upstream", 0, 500,
            cigar=[(SOFT, 1), (MATCH, 49)], seq=SEQ_FWD,
            tags=[("XM", 0), ("AS", 55)]))

        # PASS: mapped [2500, 2549) — entirely after transposon
        bam.write(aln(bam, "transp_pass_downstream", 0, 2500,
            cigar=[(SOFT, 1), (MATCH, 49)], seq=SEQ_FWD,
            tags=[("XM", 0), ("AS", 55)]))

        # PASS: mapped [2000, 2049) — starts exactly at transposon end (no overlap)
        # is_reverse=True so step 6 gets at least one reverse-strand read for rev.bw
        bam.write(aln(bam, "transp_pass_boundary", 0, 2000,
            cigar=[(MATCH, 49), (SOFT, 1)], seq=SEQ_REV,
            is_reverse=True,
            tags=[("XM", 0), ("AS", 55)]))

        # FAIL: mapped [1100, 1149) — fully inside transposon
        bam.write(aln(bam, "transp_fail_inside", 0, 1100,
            cigar=[(SOFT, 1), (MATCH, 49)], seq=SEQ_FWD,
            tags=[("XM", 0), ("AS", 55)]))

        # FAIL: mapped [970, 1019) — overlaps transposon start [1000, 1019)
        bam.write(aln(bam, "transp_fail_overlap_start", 0, 970,
            cigar=[(SOFT, 1), (MATCH, 49)], seq=SEQ_FWD,
            tags=[("XM", 0), ("AS", 55)]))

        # FAIL: mapped [1970, 2019) — overlaps transposon end [1970, 2000)
        bam.write(aln(bam, "transp_fail_overlap_end", 0, 1970,
            cigar=[(SOFT, 1), (MATCH, 49)], seq=SEQ_FWD,
            tags=[("XM", 0), ("AS", 55)]))

    write_bam("step5_input.bam", header, build)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ensure_out_dir()
    print(f"Writing test data to: {OUT_DIR}")
    print()
    print("BED files:")
    make_centromere_bed()
    make_transposon_bed()
    print()
    print("BAM files:")
    make_step1_input()
    make_step3_input()
    make_step4_input()
    make_step5_input()
    print()
    print("Done.")


if __name__ == "__main__":
    main()
