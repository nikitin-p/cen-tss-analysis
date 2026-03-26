#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  bowtie2_remap.sh -i INPUT_FASTQ_GZ -x INDEX_PREFIX -o OUTPUT_BAM [-t THREADS] [--extra-args '...']
USAGE
}

threads=20
input_fastq=""
index_prefix=""
output_bam=""
extra_args="-k 800 -D 1600 -R 1 --local"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-fastq) input_fastq="$2"; shift 2 ;;
    -x|--index-prefix) index_prefix="$2"; shift 2 ;;
    -o|--output-bam) output_bam="$2"; shift 2 ;;
    -t|--threads) threads="$2"; shift 2 ;;
    --extra-args) extra_args="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$input_fastq" ]] || { echo "Missing --input-fastq" >&2; exit 1; }
[[ -n "$index_prefix" ]] || { echo "Missing --index-prefix" >&2; exit 1; }
[[ -n "$output_bam" ]] || { echo "Missing --output-bam" >&2; exit 1; }
[[ -f "$input_fastq" ]] || { echo "Input FASTQ not found: $input_fastq" >&2; exit 1; }

mkdir -p "$(dirname "$output_bam")"

echo "[bowtie2_remap] remapping with bowtie2: $(date)"
# shellcheck disable=SC2086
bowtie2 \
  $extra_args \
  -x "$index_prefix" \
  -U "$input_fastq" \
  -p "$threads" | \
  samtools view -b - | \
  samtools sort -@ "$threads" -O bam -o "$output_bam"

samtools index -@ "$threads" "$output_bam"
echo "[bowtie2_remap] wrote $output_bam"
