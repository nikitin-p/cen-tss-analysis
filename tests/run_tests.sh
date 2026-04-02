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
        printf "  [PASS] %-45s expected %d, got %d\n" "$label" "$expected" "$got"
        ((pass++)) || true
    else
        printf "  [FAIL] %-45s expected %d, got %d\n" "$label" "$expected" "$got"
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

echo "Step 3 — pyselectal.py"
python "$SCRIPTS/pyselectal.py" \
    -n 1 -m 1 -x G \
    "$DATA/step3_input.bam" \
    "$OUT/step3_output.bam"

samtools index "$OUT/step3_output.bam"
n=$(bam_count "$OUT/step3_output.bam")
check "reads passing 5'-G filter" 3 "$n"
echo

# ── Step 4: xm_filter.sh ─────────────────────────────────────────────────────
# Input : step4_input.bam (5 reads)
# Filter: XM/read_length <= 0.10
# Expect: 2 reads pass (xm_pass_0mm, xm_pass_3mm_boundary)
#         xm_fail_4mm, xm_fail_10mm  → filtered out (>10%)
#         xm_drop_noXM               → silently dropped (no XM tag)

echo "Step 4 — xm_filter.sh"
bash "$SCRIPTS/xm_filter.sh" \
    -i "$DATA/step4_input.bam" \
    -o "$OUT/step4_output.bam" \
    -t 4

n=$(bam_count "$OUT/step4_output.bam")
check "reads passing XM filter" 2 "$n"
echo

# ── Step 5: filter_transposons.sh ─────────────────────────────────────────────
# Input : step5_input.bam (4 reads) + transposon.bed (HOR_chr1:1000-2000)
# Filter: bedtools intersect -v
# Expect: 2 reads pass (transp_pass_upstream, transp_pass_downstream)
#         transp_fail_inside, transp_fail_overlap_start → filtered out

echo "Step 5 — filter_transposons.sh"
bash "$SCRIPTS/filter_transposons.sh" \
    -i "$DATA/step5_input.bam" \
    -b "$DATA/transposon.bed" \
    -o "$OUT/step5_output.bam" \
    -t 4

n=$(bam_count "$OUT/step5_output.bam")
check "reads passing transposon filter" 2 "$n"
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
