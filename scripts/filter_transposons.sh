#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: filter_transposons.sh -i INPUT_BAM -b TRANSPOSON_BED -o OUTPUT_BAM [-t THREADS]

Removes alignments that overlap centromeric transposons using bedtools intersect -v.
Any read with any overlap to the transposon BED is discarded.
The output BAM preserves the coordinate sort order of the input.
USAGE
}

threads=20
input_bam=""
transposon_bed=""
output_bam=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-bam)    input_bam="$2";      shift 2 ;;
    -b|--bed)          transposon_bed="$2"; shift 2 ;;
    -o|--output-bam)   output_bam="$2";     shift 2 ;;
    -t|--threads)      threads="$2";        shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$input_bam" ]]      || { echo "Missing --input-bam" >&2; exit 1; }
[[ -n "$transposon_bed" ]] || { echo "Missing --bed"       >&2; exit 1; }
[[ -n "$output_bam" ]]     || { echo "Missing --output-bam" >&2; exit 1; }
[[ -f "$input_bam" ]]      || { echo "Input BAM not found: $input_bam" >&2; exit 1; }
[[ -f "$transposon_bed" ]] || { echo "Transposon BED not found: $transposon_bed" >&2; exit 1; }

mkdir -p "$(dirname "$output_bam")"

echo "[filter_transposons] removing transposon-overlapping reads: $(date)"
bedtools intersect -v -abam "$input_bam" -b "$transposon_bed" \
  | samtools view -b -o "$output_bam"

samtools index -@ "$threads" "$output_bam"
echo "[filter_transposons] wrote $output_bam"
