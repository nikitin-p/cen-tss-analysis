#!/usr/bin/env bash
# run_tests.sh — exercise pipeline steps 1, 3, 4, 5 against synthetic test data.
#
# Usage:
#   bash tests/run_tests.sh            # from repo root
#
# Prerequisites:
#   python tests/make_test_data.py     # must have been run first
#   pysam, samtools, bedtools          # must be on PATH
#
# Steps 2 (bowtie2 remap) and 6 (bamCoverage) are not tested here because they
# require an external reference index and deepTools respectively.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$REPO_ROOT/tests/data"
OUT="$REPO_ROOT/tests/out"
SCRIPTS="$REPO_ROOT/scripts"

mkdir -p "$OUT"

pass=0
fail=0

# ── helpers ───────────────────────────────────────────────────────────────────

check() {
    local label="$1"
    local expected="$2"
    local got="$3"
    if [[ "$got" -eq "$expected" ]]; then
        printf "  [PASS] %-52s expected %d, got %d\n" "$label" "$expected" "$got"
        ((pass++)) || true
    else
        printf "  [FAIL] %-52s expected %d, got %d\n" "$label" "$expected" "$got"
        ((fail++)) || true
    fi
}

# check_present BAM READNAME LABEL
#   Verifies that READNAME appears in BAM (read should have passed the filter).
check_present() {
    local bam="$1" name="$2" label="$3"
    if samtools view "$bam" | cut -f1 | grep -qx "$name"; then
        printf "  [PASS] %-52s '%s' present\n" "$label" "$name"
        ((pass++)) || true
    else
        printf "  [FAIL] %-52s '%s' missing (should be present)\n" "$label" "$name"
        ((fail++)) || true
    fi
}

# check_absent BAM READNAME LABEL
#   Verifies that READNAME does NOT appear in BAM (read should have been filtered out).
check_absent() {
    local bam="$1" name="$2" label="$3"
    if ! samtools view "$bam" | cut -f1 | grep -qx "$name"; then
        printf "  [PASS] %-52s '%s' absent\n" "$label" "$name"
        ((pass++)) || true
    else
        printf "  [FAIL] %-52s '%s' present (should be absent)\n" "$label" "$name"
        ((fail++)) || true
    fi
}

bam_count() {
    samtools view -c "$1" 2>/dev/null
}

fastq_count() {
    # count reads: divide line count by 4
    zcat "$1" 2>/dev/null | wc -l | awk '{print $1/4}'
}

# ── preflight: verify test data exists ────────────────────────────────────────

required_files=(
    "$DATA/centromere.bed"
    "$DATA/transposon.bed"
    "$DATA/step1_input.bam"
    "$DATA/step3_input.bam"
    "$DATA/step4_input.bam"
    "$DATA/step5_input.bam"
)

echo "Checking test data..."
for f in "${required_files[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "  MISSING: $f"
        echo "  Run:  python tests/make_test_data.py"
        exit 1
    fi
done
echo "  All test data present."
echo

# ── Step 1: get_remap.sh ──────────────────────────────────────────────────────
# Input : step1_input.bam + centromere.bed
# Expect: 5 centromeric reads + 3 unmapped = 8 reads in output FASTQ.

echo "Step 1 — get_remap.sh"
bash "$SCRIPTS/get_remap.sh" \
    -i "$DATA/step1_input.bam" \
    -b "$DATA/centromere.bed" \
    -o "$OUT/step1_output.fastq.gz" \
    -t 4

n=$(fastq_count "$OUT/step1_output.fastq.gz")
check "reads in FASTQ (cen+unmapped)" 8 "$n"
echo

# ── Step 3: pyselectal.py ─────────────────────────────────────────────────────
# Input : step3_input.bam (7 reads)
# Filter: -n 1 -m 1 -x G  (exactly 1 5'-soft-clipped G)
# Expect: 3 reads pass (fwd_pass_1G, fwd_pass_2G, rev_pass_1C)
#         4 reads fail (wrong base, wrong clip length, no clip)

echo "Step 3 — pyselectal.py"
python "$SCRIPTS/pyselectal.py" \
    -n 1 -m 1 -x G \
    "$DATA/step3_input.bam" \
    "$OUT/step3_output.bam"

samtools index "$OUT/step3_output.bam"
n=$(bam_count "$OUT/step3_output.bam")
check "total reads passing 5'-G filter" 3 "$n"

echo "  -- per-read verification:"
BAM3="$OUT/step3_output.bam"
check_present "$BAM3" "fwd_pass_1G"      "fwd 1S+50M, first base G"
check_present "$BAM3" "fwd_pass_2G"      "fwd 1S+50M, first base G (2nd)"
check_present "$BAM3" "rev_pass_1C"      "rev 50M+1S, last base C (revcomp G)"
check_absent  "$BAM3" "fwd_fail_notG"    "fwd 1S+50M but first base A"
check_absent  "$BAM3" "fwd_fail_2Sclip"  "fwd 2S+49M (2-bp clip, not 1)"
check_absent  "$BAM3" "fwd_fail_noSclip" "fwd 51M (no soft-clip)"
check_absent  "$BAM3" "rev_fail_notC"    "rev 50M+1S but last base A, not C"
echo

# ── Step 4: xm_filter.sh ─────────────────────────────────────────────────────
# Input : step4_input.bam (5 reads, all 30 bp)
# Filter: XM/read_length <= 0.10  (i.e. >= 90% match rate)
# Expect: 2 reads pass
#   xm_pass_0mm           XM=0,  0/30=0.000 — no mismatches
#   xm_pass_3mm_boundary  XM=3,  3/30=0.100 — exactly at the 10% boundary
# 3 reads filtered/dropped:
#   xm_fail_4mm           XM=4,  4/30=0.133 — just over 10%
#   xm_fail_10mm          XM=10, 10/30=0.333 — far over 10%
#   xm_drop_noXM          no XM tag — silently dropped by awk

echo "Step 4 — xm_filter.sh"
bash "$SCRIPTS/xm_filter.sh" \
    -i "$DATA/step4_input.bam" \
    -o "$OUT/step4_output.bam" \
    -t 4

n=$(bam_count "$OUT/step4_output.bam")
check "total reads passing XM filter" 2 "$n"

echo "  -- per-read verification:"
BAM4="$OUT/step4_output.bam"
check_present "$BAM4" "xm_pass_0mm"           "XM=0  (0/30=0.000)"
check_present "$BAM4" "xm_pass_3mm_boundary"  "XM=3  (3/30=0.100, at boundary)"
check_absent  "$BAM4" "xm_fail_4mm"           "XM=4  (4/30=0.133, >10%)"
check_absent  "$BAM4" "xm_fail_10mm"          "XM=10 (10/30=0.333, >10%)"
check_absent  "$BAM4" "xm_drop_noXM"          "no XM tag (silently dropped)"
echo

# ── Step 5: filter_transposons.sh ─────────────────────────────────────────────
# Input : step5_input.bam (6 reads) + transposon.bed (HOR_chr1 [1000,2000))
# Filter: bedtools intersect -v  (removes any read overlapping transposon BED)
# Expect: 3 reads pass
#   transp_pass_upstream   mapped [500,  549) — before transposon
#   transp_pass_downstream mapped [2500, 2549) — after transposon
#   transp_pass_boundary   mapped [2000, 2049) — starts exactly at transposon end
# 3 reads filtered:
#   transp_fail_inside        mapped [1100, 1149) — fully inside
#   transp_fail_overlap_start mapped [970,  1019) — overlaps transposon start
#   transp_fail_overlap_end   mapped [1970, 2019) — overlaps transposon end

echo "Step 5 — filter_transposons.sh"
bash "$SCRIPTS/filter_transposons.sh" \
    -i "$DATA/step5_input.bam" \
    -b "$DATA/transposon.bed" \
    -o "$OUT/step5_output.bam" \
    -t 4

n=$(bam_count "$OUT/step5_output.bam")
check "total reads passing transposon filter" 3 "$n"

echo "  -- per-read verification:"
BAM5="$OUT/step5_output.bam"
check_present "$BAM5" "transp_pass_upstream"      "mapped [500,549)  before transposon"
check_present "$BAM5" "transp_pass_downstream"    "mapped [2500,2549) after transposon"
check_present "$BAM5" "transp_pass_boundary"      "mapped [2000,2049) starts at TE end"
check_absent  "$BAM5" "transp_fail_inside"        "mapped [1100,1149) fully inside TE"
check_absent  "$BAM5" "transp_fail_overlap_start" "mapped [970,1019)  overlaps TE start"
check_absent  "$BAM5" "transp_fail_overlap_end"   "mapped [1970,2019) overlaps TE end"
echo

# ── Step 6: bamcoverage_5p.sh ─────────────────────────────────────────────────
# Input : step5_output.bam (3 reads: 2 forward, 1 reverse)
# Expect: step6_output.fwd.bw and step6_output.rev.bw created and non-empty

echo "Step 6 — bamcoverage_5p.sh"
bash "$SCRIPTS/bamcoverage_5p.sh" \
    -i "$OUT/step5_output.bam" \
    -o "$OUT/step6_output" \
    -t 4

check_file() {
    local label="$1"
    local path="$2"
    if [[ -f "$path" && -s "$path" ]]; then
        printf "  [PASS] %-45s exists and non-empty\n" "$label"
        ((pass++)) || true
    else
        printf "  [FAIL] %-45s missing or empty\n" "$label"
        ((fail++)) || true
    fi
}

check_file "step6_output.fwd.bw" "$OUT/step6_output.fwd.bw"
check_file "step6_output.rev.bw" "$OUT/step6_output.rev.bw"
echo

# ── summary ───────────────────────────────────────────────────────────────────

total=$((pass + fail))
echo "────────────────────────────────────────────────────────"
echo "Results: $pass / $total passed"
if [[ "$fail" -gt 0 ]]; then
    echo "FAILED: $fail test(s)"
    exit 1
else
    echo "All tests passed."
fi
