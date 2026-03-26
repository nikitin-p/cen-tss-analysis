#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  xm_filter.sh -i INPUT_BAM -o OUTPUT_BAM [-t THREADS] [--min-match-frac 0.90]
USAGE
}

threads=20
input_bam=""
output_bam=""
min_match_frac="0.90"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-bam) input_bam="$2"; shift 2 ;;
    -o|--output-bam) output_bam="$2"; shift 2 ;;
    -t|--threads) threads="$2"; shift 2 ;;
    --min-match-frac) min_match_frac="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$input_bam" ]] || { echo "Missing --input-bam" >&2; exit 1; }
[[ -n "$output_bam" ]] || { echo "Missing --output-bam" >&2; exit 1; }
[[ -f "$input_bam" ]] || { echo "Input BAM not found: $input_bam" >&2; exit 1; }

mkdir -p "$(dirname "$output_bam")"

awk_threshold=$(awk -v x="$min_match_frac" 'BEGIN{printf "%.12f", 1-x}')

echo "[xm_filter] filtering by XM/read_length <= $awk_threshold : $(date)"
samtools view -h "$input_bam" | \
  awk -v max_mm_frac="$awk_threshold" 'BEGIN{OFS="\t"}
    /^@/ {print; next}
    {
      readlen = length($10)
      if (readlen == 0) next
      for (i = 12; i <= NF; i++) {
        if ($i ~ /^XM:i:/) {
          split($i, xm, ":")
          if ((xm[3] / readlen) <= max_mm_frac) print
          break
        }
      }
    }' | \
  samtools view -b -o "$output_bam"

samtools index -@ "$threads" "$output_bam"
echo "[xm_filter] wrote $output_bam"
