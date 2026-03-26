#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<USAGE
Usage: select_tss_alignments.sh -i INPUT_BAM -o OUTPUT_BAM [-t THREADS]
                                 [--pyselectal-args '...']

Selects alignments with a 5'-end G soft clip using pyselectal.py.

Default filter (-n 1 -m 1 -x G):
  Keep reads with exactly one 5'-soft-clipped base that is G (forward strand)
  or C (reverse strand). Both the CIGAR pattern (1S prefix) and the actual
  nucleotide are checked. This is the biologically correct criterion for CAGE
  TSS identification.

Note on old behaviour:
  The original get5GsoftclipBAM.sh only checked the CIGAR pattern (1S prefix)
  without verifying the base identity. This script checks both, which is more
  precise. The difference affects reads where the soft-clipped base is not G.
USAGE
}

threads=1
input_bam=""
output_bam=""
pyselectal_args="-n 1 -m 1 -x G"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-bam)       input_bam="$2";       shift 2 ;;
    -o|--output-bam)      output_bam="$2";      shift 2 ;;
    -t|--threads)         threads="$2";         shift 2 ;;
    --pyselectal-args)    pyselectal_args="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$input_bam" ]]  || { echo "Missing --input-bam"  >&2; exit 1; }
[[ -n "$output_bam" ]] || { echo "Missing --output-bam" >&2; exit 1; }
[[ -f "$input_bam" ]]  || { echo "Input BAM not found: $input_bam" >&2; exit 1; }

mkdir -p "$(dirname "$output_bam")"

echo "[select_tss_alignments] selecting 5'-G alignments: $(date)"
# shellcheck disable=SC2086
python3 "${SCRIPT_DIR}/pyselectal.py" \
  $pyselectal_args \
  -t "$threads" \
  "$input_bam" \
  "$output_bam"

samtools index -@ "$threads" "$output_bam"
echo "[select_tss_alignments] wrote $output_bam"
