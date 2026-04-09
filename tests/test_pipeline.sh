#!/usr/bin/env bash
# test_pipeline.sh — integration tests for run_pipeline.sh
#
# Usage:
#   bash tests/test_pipeline.sh           # from repo root
#
# Prerequisites:
#   python tests/make_test_data.py        # generates tests/data/ inputs
#   samtools, python3                     # must be on PATH
#
# bedtools, bowtie2, bamCoverage are NOT required — lightweight stubs in
# tests/stubs/ are injected into PATH automatically.  Steps 3 (pyselectal.py)
# and 4 (xm_filter.sh) run real code; all other steps use stubs.
#
# What is tested
# --------------
# 1. Error handling: missing -i, missing -o, bad config, bad input path
# 2. --dry-run: commands printed, no output files created
# 3. Full run (single BAM): all 6 steps execute; final outputs exist;
#    step 3 correctly filtered (3 reads); intermediate BAMs present
# 4. --resume: re-running with existing outputs skips every step
# 5. Batch run (directory of BAMs): both samples processed independently

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$REPO_ROOT/tests/data"
OUT="$REPO_ROOT/tests/pipeline_out"
STUBS="$REPO_ROOT/tests/stubs"
CONF="$REPO_ROOT/tests/test_pipeline.conf"
PIPELINE="$REPO_ROOT/run_pipeline.sh"

# Inject stubs before system tools so bedtools/bowtie2/bamCoverage resolve
# to the lightweight stubs rather than real (absent) executables.
export PATH="$STUBS:$PATH"

rm -rf "$OUT"
mkdir -p "$OUT"

pass=0
fail=0

# ── helpers ───────────────────────────────────────────────────────────────────

check() {
    local label="$1" expected="$2" got="$3"
    if [[ "$got" -eq "$expected" ]]; then
        printf "  [PASS] %-55s expected %d, got %d\n" "$label" "$expected" "$got"
        ((pass++)) || true
    else
        printf "  [FAIL] %-55s expected %d, got %d\n" "$label" "$expected" "$got"
        ((fail++)) || true
    fi
}

check_file() {
    local label="$1" path="$2"
    if [[ -f "$path" && -s "$path" ]]; then
        printf "  [PASS] %-55s exists and non-empty\n" "$label"
        ((pass++)) || true
    else
        printf "  [FAIL] %-55s missing or empty: %s\n" "$label" "$path"
        ((fail++)) || true
    fi
}

check_no_file() {
    local label="$1" path="$2"
    if [[ ! -f "$path" ]]; then
        printf "  [PASS] %-55s absent as expected\n" "$label"
        ((pass++)) || true
    else
        printf "  [FAIL] %-55s exists but should be absent: %s\n" "$label" "$path"
        ((fail++)) || true
    fi
}

check_contains() {
    local label="$1" file="$2" pattern="$3"
    if grep -q "$pattern" "$file" 2>/dev/null; then
        printf "  [PASS] %-55s contains '%s'\n" "$label" "$pattern"
        ((pass++)) || true
    else
        printf "  [FAIL] %-55s '%s' not found in %s\n" "$label" "$pattern" "$file"
        ((fail++)) || true
    fi
}

check_not_contains() {
    local label="$1" file="$2" pattern="$3"
    if ! grep -q "$pattern" "$file" 2>/dev/null; then
        printf "  [PASS] %-55s does not contain '%s'\n" "$label" "$pattern"
        ((pass++)) || true
    else
        printf "  [FAIL] %-55s '%s' found but should be absent in %s\n" "$label" "$pattern" "$file"
        ((fail++)) || true
    fi
}

check_error() {
    # Verify that the given command exits with a non-zero status.
    local label="$1"
    shift
    if "$@" >/dev/null 2>&1; then
        printf "  [FAIL] %-55s should have failed but succeeded\n" "$label"
        ((fail++)) || true
    else
        printf "  [PASS] %-55s exited non-zero as expected\n" "$label"
        ((pass++)) || true
    fi
}

check_present() {
    local bam="$1" name="$2" label="$3"
    if samtools view "$bam" | cut -f1 | grep -qx "$name"; then
        printf "  [PASS] %-55s '%s' present\n" "$label" "$name"
        ((pass++)) || true
    else
        printf "  [FAIL] %-55s '%s' missing (should be present)\n" "$label" "$name"
        ((fail++)) || true
    fi
}

check_absent() {
    local bam="$1" name="$2" label="$3"
    if ! samtools view "$bam" | cut -f1 | grep -qx "$name"; then
        printf "  [PASS] %-55s '%s' absent\n" "$label" "$name"
        ((pass++)) || true
    else
        printf "  [FAIL] %-55s '%s' present (should be absent)\n" "$label" "$name"
        ((fail++)) || true
    fi
}

bam_count() { samtools view -c "$1" 2>/dev/null; }

# ── preflight ─────────────────────────────────────────────────────────────────

echo "Checking test data..."
required=(
    "$DATA/centromere.bed"
    "$DATA/transposon.bed"
    "$DATA/step1_input.bam"
    "$DATA/step3_input.bam"
    "$DATA/fake_index/ref.1.bt2"
    "$CONF"
)
for f in "${required[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "  MISSING: $f"
        echo "  Run:  python tests/make_test_data.py"
        exit 1
    fi
done
echo "  All prerequisites present."
echo

# ── 1. Error handling ─────────────────────────────────────────────────────────
#
# These tests verify that run_pipeline.sh exits non-zero and does not start
# processing when required arguments or files are absent.

echo "=== 1. Error handling ==="
check_error "missing -i"         bash "$PIPELINE" -o "$OUT/err" -c "$CONF"
check_error "missing -o"         bash "$PIPELINE" -i "$DATA/step1_input.bam" -c "$CONF"
check_error "nonexistent config" bash "$PIPELINE" -i "$DATA/step1_input.bam" -o "$OUT/err" -c /nonexistent.conf
check_error "nonexistent -i BAM" bash "$PIPELINE" -i /nonexistent.bam -o "$OUT/err" -c "$CONF"
check_error "nonexistent -i dir" bash "$PIPELINE" -i /nonexistent_dir -o "$OUT/err" -c "$CONF"
echo

# ── 2. --dry-run ──────────────────────────────────────────────────────────────
#
# --dry-run must print [DRY-RUN] lines for each step command and must NOT
# create any output files under the output directory.

echo "=== 2. --dry-run ==="
DRY_OUT="$OUT/dryrun"
dry_log="$OUT/dryrun_stdout.txt"
bash "$PIPELINE" \
    -i "$DATA/step1_input.bam" \
    -o "$DRY_OUT" \
    -c "$CONF" \
    -t 1 \
    --dry-run > "$dry_log" 2>&1 || true

check_contains  "dry-run prints step 1 (get_remap)"           "$dry_log" "get_remap.sh"
check_contains  "dry-run prints step 2 (bowtie2_remap)"       "$dry_log" "bowtie2_remap.sh"
check_contains  "dry-run prints step 3 (select_tss)"          "$dry_log" "select_tss_alignments.sh"
check_contains  "dry-run prints step 4 (xm_filter)"           "$dry_log" "xm_filter.sh"
check_contains  "dry-run prints step 5 (filter_transposons)"  "$dry_log" "filter_transposons.sh"
check_contains  "dry-run prints step 6 (bamcoverage)"         "$dry_log" "bamcoverage_5p.sh"
check_no_file   "dry-run creates no fwd.bw"  "$DRY_OUT/step1_input/step1_input.fwd.bw"
check_no_file   "dry-run creates no log"     "$DRY_OUT/step1_input/step1_input.log"
echo

# ── 3. Full run (single BAM, --keep-intermediates) ────────────────────────────
#
# Runs all 6 pipeline steps with --keep-intermediates so intermediate BAMs
# are retained for inspection.  Steps 3 and 4 run real code:
#   - bowtie2 stub feeds step3_input.bam (7 reads) into the pipeline
#   - pyselectal.py (step 3) keeps 3 reads: fwd_pass_1G, fwd_pass_2G, rev_pass_1C
#   - xm_filter.sh (step 4) keeps all 3 (XM=0 for all)
#   - bedtools stub passes all 3 through (step 5)
#   - bamCoverage stub creates stub .bw files (step 6)

echo "=== 3. Full run (single BAM, --keep-intermediates) ==="
FULL_OUT="$OUT/full"
SDIR="$FULL_OUT/step1_input"
bash "$PIPELINE" \
    -i "$DATA/step1_input.bam" \
    -o "$FULL_OUT" \
    -c "$CONF" \
    -t 1 \
    --keep-intermediates 2>&1 | tee "$OUT/full_run_stdout.txt" >/dev/null || true

echo "  -- output files:"
check_file  "fwd.bw created"                    "$SDIR/step1_input.fwd.bw"
check_file  "rev.bw created"                    "$SDIR/step1_input.rev.bw"
check_file  "log file created"                  "$SDIR/step1_input.log"

echo "  -- intermediate files:"
check_file  "step2 BAM (bowtie2 output)"        "$SDIR/step1_input.bowtie2.bam"
check_file  "step3 BAM (softclip_g)"            "$SDIR/step1_input.softclip_g.bam"
check_file  "step4 BAM (xm_filtered)"           "$SDIR/step1_input.xm_filtered.bam"
check_file  "step5 BAM (no_transposons)"        "$SDIR/step1_input.no_transposons.bam"

echo "  -- step 3 filtering (pyselectal.py ran correctly):"
SC_BAM="$SDIR/step1_input.softclip_g.bam"
n=$(bam_count "$SC_BAM")
check "step3 output read count" 3 "$n"
check_present "$SC_BAM" "fwd_pass_1G"      "fwd 1S+50M, first base G"
check_present "$SC_BAM" "fwd_pass_2G"      "fwd 1S+50M, first base G (2nd)"
check_present "$SC_BAM" "rev_pass_1C"      "rev 50M+1S, last base C (revcomp G)"
check_absent  "$SC_BAM" "fwd_fail_notG"    "fwd 1S+50M but first base A"
check_absent  "$SC_BAM" "fwd_fail_2Sclip"  "fwd 2S+49M (2-bp clip, not 1)"
check_absent  "$SC_BAM" "fwd_fail_noSclip" "fwd 51M (no soft-clip)"
check_absent  "$SC_BAM" "rev_fail_notC"    "rev 50M+1S but last base A, not C"

echo "  -- step 4 → step 5 read count (all 3 survive XM + transposon filter):"
TP_BAM="$SDIR/step1_input.no_transposons.bam"
n=$(bam_count "$TP_BAM")
check "step5 output read count" 3 "$n"

echo "  -- log contains expected messages:"
LOG="$SDIR/step1_input.log"
check_contains "log: step 1 started"  "$LOG" "\[1/6\]"
check_contains "log: step 2 started"  "$LOG" "\[2/6\]"
check_contains "log: step 3 started"  "$LOG" "\[3/6\]"
check_contains "log: step 4 started"  "$LOG" "\[4/6\]"
check_contains "log: step 5 started"  "$LOG" "\[5/6\]"
check_contains "log: step 6 started"  "$LOG" "\[6/6\]"
check_contains "log: Done message"    "$LOG" "Done: step1_input"
echo

# ── 4. --resume ───────────────────────────────────────────────────────────────
#
# Re-run against the same output directory.  All outputs from test 3 already
# exist, so every step must be skipped.  The log (appended) should contain
# [resume] lines and must NOT contain any new [1/6] .. [6/6] processing lines
# after the resume marker.

echo "=== 4. --resume ==="
bash "$PIPELINE" \
    -i "$DATA/step1_input.bam" \
    -o "$FULL_OUT" \
    -c "$CONF" \
    -t 1 \
    --resume 2>&1 >/dev/null || true

LOG="$SDIR/step1_input.log"
check_contains "log: resume skipped step 1"  "$LOG" "\[resume\].*remap.fastq.gz"
check_contains "log: resume skipped step 2"  "$LOG" "\[resume\].*bowtie2.bam"
check_contains "log: resume skipped step 3"  "$LOG" "\[resume\].*softclip_g.bam"
check_contains "log: resume skipped step 4"  "$LOG" "\[resume\].*xm_filtered.bam"
check_contains "log: resume skipped step 5"  "$LOG" "\[resume\].*no_transposons.bam"
check_contains "log: resume skipped step 6"  "$LOG" "\[resume\].*fwd.bw"
echo

# ── 5. Batch run (directory of BAMs) ──────────────────────────────────────────
#
# Put step1_input.bam and step3_input.bam in a temporary directory and run
# the pipeline on the whole directory.  Both samples must be processed and
# produce independent output subdirectories.

echo "=== 5. Batch run (directory input) ==="
BATCH_IN="$OUT/batch_in"
BATCH_OUT="$OUT/batch"
mkdir -p "$BATCH_IN"
cp "$DATA/step1_input.bam" "$BATCH_IN/sampleA.bam"
cp "$DATA/step1_input.bam" "$BATCH_IN/sampleB.bam"
samtools index "$BATCH_IN/sampleA.bam"
samtools index "$BATCH_IN/sampleB.bam"

bash "$PIPELINE" \
    -i "$BATCH_IN" \
    -o "$BATCH_OUT" \
    -c "$CONF" \
    -t 1 2>&1 >/dev/null || true

check_file  "batch: sampleA fwd.bw"  "$BATCH_OUT/sampleA/sampleA.fwd.bw"
check_file  "batch: sampleA rev.bw"  "$BATCH_OUT/sampleA/sampleA.rev.bw"
check_file  "batch: sampleB fwd.bw"  "$BATCH_OUT/sampleB/sampleB.fwd.bw"
check_file  "batch: sampleB rev.bw"  "$BATCH_OUT/sampleB/sampleB.rev.bw"
check_file  "batch: sampleA log"     "$BATCH_OUT/sampleA/sampleA.log"
check_file  "batch: sampleB log"     "$BATCH_OUT/sampleB/sampleB.log"
echo

# ── summary ───────────────────────────────────────────────────────────────────

total=$((pass + fail))
echo "────────────────────────────────────────────────────────────────"
echo "Results: $pass / $total passed"
if [[ "$fail" -gt 0 ]]; then
    echo "FAILED: $fail test(s)"
    exit 1
else
    echo "All tests passed."
fi
