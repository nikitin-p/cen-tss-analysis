#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CONFIG_FILE="${SCRIPT_DIR}/pipeline.conf"

usage() {
  cat <<USAGE
Usage: run_pipeline.sh -i INPUT -o OUTPUT_DIR [options]

Required:
  -i, --input PATH         BAM file or directory containing BAM files
  -o, --output DIR         Root output directory

Options:
  -c, --config FILE        Config file (default: pipeline.conf next to this script)
  -t, --threads N          Thread count (overrides THREADS in config)
  --keep-intermediates     Keep intermediate BAM and FASTQ files (default: delete)
  --resume                 Skip steps whose output already exists
  --dry-run                Print commands without executing anything
  -h, --help               Show this message

Output per sample (under OUTPUT_DIR/SAMPLE/):
  SAMPLE.fwd.bw    5-prime coverage bigWig, forward strand (+)
  SAMPLE.rev.bw    5-prime coverage bigWig, reverse strand (-)
  SAMPLE.log       Full run log for this sample

Example:
  run_pipeline.sh -i data/bams/ -o results/ -c pipeline.conf -t 16
  run_pipeline.sh -i data/sampleA.bam -o results/ --resume
USAGE
}

INPUT=""
OUTDIR=""
THREADS_OVERRIDE=""
KEEP_INTERMEDIATES=0
RESUME=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)            INPUT="$2";             shift 2 ;;
    -o|--output)           OUTDIR="$2";            shift 2 ;;
    -c|--config)           CONFIG_FILE="$2";       shift 2 ;;
    -t|--threads)          THREADS_OVERRIDE="$2";  shift 2 ;;
    --keep-intermediates)  KEEP_INTERMEDIATES=1;   shift   ;;
    --resume)              RESUME=1;               shift   ;;
    --dry-run)             DRY_RUN=1;              shift   ;;
    -h|--help)             usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$INPUT" ]]  || { echo "ERROR: -i/--input is required"  >&2; usage; exit 1; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: -o/--output is required" >&2; usage; exit 1; }
[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: config not found: $CONFIG_FILE" >&2; exit 1; }

# shellcheck disable=SC1090
source "$CONFIG_FILE"

[[ -n "${THREADS_OVERRIDE:-}" ]] && THREADS="$THREADS_OVERRIDE"

# ── helpers ───────────────────────────────────────────────────────────────────

SAMPLE_LOG=""   # set to per-sample log path inside process_sample

log() {
  local msg="[$(date '+%F %T')] $*"
  echo "$msg"
  [[ -n "${SAMPLE_LOG:-}" ]] && echo "$msg" >> "$SAMPLE_LOG"
}

err() { echo "ERROR: $*" >&2; exit 1; }

need() { command -v "$1" >/dev/null 2>&1 || err "Required tool not found: $1"; }

run_cmd() {
  # Dry-run: print command to screen.
  # Normal:  execute and append stdout+stderr to the per-sample log.
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '[DRY-RUN]'
    printf ' %q' "$@"
    printf '\n'
    return 0
  fi
  if [[ -n "${SAMPLE_LOG:-}" ]]; then
    "$@" >> "$SAMPLE_LOG" 2>&1
  else
    "$@"
  fi
}

validate_output() {
  [[ "$DRY_RUN" -eq 1 ]] && return 0
  [[ -s "$1" ]] || err "Step output was not created or is empty: $1"
}

safe_rm() {
  # Delete files unless --keep-intermediates is set.
  [[ "$KEEP_INTERMEDIATES" -eq 1 ]] && return 0
  local f
  for f in "$@"; do
    if [[ "$DRY_RUN" -eq 1 ]]; then
      echo "[DRY-RUN] rm -f $f"
    elif [[ -f "$f" ]]; then
      rm -f "$f"
    fi
  done
}

step_is_done() {
  # Returns 0 (skip) when --resume is active and the output already exists.
  if [[ "$RESUME" -eq 1 && -s "$1" ]]; then
    log "  [resume] skipping, output exists: $(basename "$1")"
    return 0
  fi
  return 1
}

# ── tool validation ───────────────────────────────────────────────────────────

need samtools
need bedtools
need bowtie2
need bamCoverage
need python3

# ── script validation ─────────────────────────────────────────────────────────

for _s in get_remap.sh bowtie2_remap.sh select_tss_alignments.sh \
          xm_filter.sh filter_transposons.sh bamcoverage_5p.sh; do
  [[ -f "${SCRIPT_DIR}/scripts/${_s}" ]] || err "Missing script: scripts/${_s}"
done
[[ -f "$PYSELECTAL_PY" ]] || err "pyselectal.py not found: $PYSELECTAL_PY"

# ── config validation ─────────────────────────────────────────────────────────

[[ -f "$CENTROMERE_BED" ]] || err "CENTROMERE_BED not found: $CENTROMERE_BED"
[[ -f "$TRANSPOSON_BED" ]] || err "TRANSPOSON_BED not found: $TRANSPOSON_BED"
[[ -f "${BOWTIE2_INDEX_PREFIX}.1.bt2" ]] \
  || [[ -f "${BOWTIE2_INDEX_PREFIX}.1.bt2l" ]] \
  || err "Bowtie2 index not found at prefix: $BOWTIE2_INDEX_PREFIX"

# ── collect input BAMs ────────────────────────────────────────────────────────

declare -a BAM_LIST=()

if [[ -f "$INPUT" ]]; then
  [[ "$INPUT" == *.bam ]] || err "Input file does not look like a BAM: $INPUT"
  BAM_LIST=("$INPUT")
elif [[ -d "$INPUT" ]]; then
  while IFS= read -r -d '' bam; do
    BAM_LIST+=("$bam")
  done < <(find "$INPUT" -maxdepth 1 -name '*.bam' -print0 | sort -z)
  [[ ${#BAM_LIST[@]} -gt 0 ]] || err "No .bam files found in: $INPUT"
else
  err "Input not found (not a file or directory): $INPUT"
fi

mkdir -p "$OUTDIR"

# ── per-sample processing ─────────────────────────────────────────────────────

process_sample() {
  local input_bam="$1"
  local sample="$2"
  local sdir="$OUTDIR/$sample"
  SAMPLE_LOG="$sdir/${sample}.log"
  mkdir -p "$sdir"

  log "=== Sample: $sample ==="
  log "    input:  $input_bam"
  log "    output: $sdir"

  local fastq="${sdir}/${sample}.remap.fastq.gz"
  local bowtie_bam="${sdir}/${sample}.bowtie2.bam"
  local softclip_bam="${sdir}/${sample}.softclip_g.bam"
  local xm_bam="${sdir}/${sample}.xm_filtered.bam"
  local no_tp_bam="${sdir}/${sample}.no_transposons.bam"
  local bw_fwd="${sdir}/${sample}.fwd.bw"
  local bw_rev="${sdir}/${sample}.rev.bw"

  # Step 1 ── extract centromeric and unmapped reads → FASTQ for remapping
  if ! step_is_done "$fastq"; then
    log "  [1/6] extracting centromeric and unmapped reads"
    run_cmd bash "${SCRIPT_DIR}/scripts/get_remap.sh" \
      -i "$input_bam" \
      -b "$CENTROMERE_BED" \
      -o "$fastq" \
      -t "$THREADS"
    validate_output "$fastq"
  fi

  # Step 2 ── remap with bowtie2 (multimappers retained)
  if ! step_is_done "$bowtie_bam"; then
    log "  [2/6] remapping with bowtie2"
    run_cmd bash "${SCRIPT_DIR}/scripts/bowtie2_remap.sh" \
      -i "$fastq" \
      -x "$BOWTIE2_INDEX_PREFIX" \
      -o "$bowtie_bam" \
      -t "$THREADS" \
      --extra-args "${BOWTIE2_EXTRA_ARGS[*]}"
    validate_output "$bowtie_bam"
  fi
  safe_rm "$fastq"

  # Step 3 ── select alignments with exactly one 5'-soft-clipped G
  if ! step_is_done "$softclip_bam"; then
    log "  [3/6] selecting 5'-G soft-clipped alignments"
    run_cmd bash "${SCRIPT_DIR}/scripts/select_tss_alignments.sh" \
      -i "$bowtie_bam" \
      -o "$softclip_bam" \
      -t "$THREADS" \
      --pyselectal-args "${PYSELECTAL_ARGS[*]}"
    validate_output "$softclip_bam"
  fi
  safe_rm "$bowtie_bam" "${bowtie_bam}.bai"

  # Step 4 ── filter by mismatch rate (XM tag < 10% mismatches)
  if ! step_is_done "$xm_bam"; then
    log "  [4/6] filtering by mismatch rate (XM)"
    run_cmd bash "${SCRIPT_DIR}/scripts/xm_filter.sh" \
      -i "$softclip_bam" \
      -o "$xm_bam" \
      -t "$THREADS" \
      --min-match-frac "$MIN_MATCH_FRAC"
    validate_output "$xm_bam"
  fi
  safe_rm "$softclip_bam" "${softclip_bam}.bai"

  # Step 5 ── remove reads overlapping centromeric transposons
  if ! step_is_done "$no_tp_bam"; then
    log "  [5/6] filtering transposon-overlapping alignments"
    run_cmd bash "${SCRIPT_DIR}/scripts/filter_transposons.sh" \
      -i "$xm_bam" \
      -b "$TRANSPOSON_BED" \
      -o "$no_tp_bam" \
      -t "$THREADS"
    validate_output "$no_tp_bam"
  fi
  safe_rm "$xm_bam" "${xm_bam}.bai"

  # Step 6 ── generate strand-specific 5'-end bigWigs
  if ! step_is_done "$bw_fwd"; then
    log "  [6/6] generating strand-specific 5'-end bigWigs"
    run_cmd bash "${SCRIPT_DIR}/scripts/bamcoverage_5p.sh" \
      -i "$no_tp_bam" \
      -o "${sdir}/${sample}" \
      -t "$THREADS"
    validate_output "$bw_fwd"
    validate_output "$bw_rev"
  fi
  # Final BAM is deleted by default; keep it only with --keep-intermediates.
  safe_rm "$no_tp_bam" "${no_tp_bam}.bai"

  log "  => $(basename "$bw_fwd")"
  log "  => $(basename "$bw_rev")"
  log "=== Done: $sample ==="
}

# ── main loop ─────────────────────────────────────────────────────────────────

declare -a FAILED_SAMPLES=()

for _bam in "${BAM_LIST[@]}"; do
  _sample=$(basename "$_bam" .bam)
  # Run process_sample in a subshell so a failure in one sample does not abort
  # the entire batch; set -e errors inside the subshell exit only the subshell.
  if ( process_sample "$_bam" "$_sample" ); then
    log "Finished: $_sample"
  else
    log "FAILED:   $_sample  (see $OUTDIR/$_sample/$_sample.log)"
    FAILED_SAMPLES+=("$_sample")
  fi
done

if [[ ${#FAILED_SAMPLES[@]} -gt 0 ]]; then
  echo "" >&2
  echo "ERROR: the following samples failed: ${FAILED_SAMPLES[*]}" >&2
  exit 1
fi

log "All samples complete."
