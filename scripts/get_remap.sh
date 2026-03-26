#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  get_remap.sh -i INPUT_BAM -b CENTROMERE_BED -o OUTPUT_FASTQ_GZ [-t THREADS] [--keep-intermediates]
USAGE
}

threads=20
keep_intermediates=0
input_bam=""
centromere_bed=""
output_fastq_gz=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-bam) input_bam="$2"; shift 2 ;;
    -b|--bed) centromere_bed="$2"; shift 2 ;;
    -o|--output-fastq) output_fastq_gz="$2"; shift 2 ;;
    -t|--threads) threads="$2"; shift 2 ;;
    --keep-intermediates) keep_intermediates=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$input_bam" ]] || { echo "Missing --input-bam" >&2; exit 1; }
[[ -n "$centromere_bed" ]] || { echo "Missing --bed" >&2; exit 1; }
[[ -n "$output_fastq_gz" ]] || { echo "Missing --output-fastq" >&2; exit 1; }
[[ -f "$input_bam" ]] || { echo "Input BAM not found: $input_bam" >&2; exit 1; }
[[ -f "$centromere_bed" ]] || { echo "Centromere BED not found: $centromere_bed" >&2; exit 1; }

outdir=$(dirname "$output_fastq_gz")
mkdir -p "$outdir"

prefix=$(basename "$output_fastq_gz")
prefix=${prefix%.gz}
prefix=${prefix%.fastq}
prefix=${prefix%.fq}
workdir=$(mktemp -d "${TMPDIR:-/tmp}/get_remap.XXXXXX")
trap 'rm -rf "$workdir"' EXIT

unmapped_bam="$workdir/${prefix}.unmapped.bam"
mapped_bam="$workdir/${prefix}.mapped.bam"
centromere_bam="$workdir/${prefix}.centromere.bam"
unmapped_fastq="$workdir/${prefix}.unmapped.fastq"
centromere_fastq="$workdir/${prefix}.centromere.fastq"
merged_fastq="$workdir/${prefix}.remap.fastq"

echo "[get_remap] extracting unmapped reads: $(date)"
samtools view -@ "$threads" -b -f 4 "$input_bam" > "$unmapped_bam"

echo "[get_remap] extracting mapped reads: $(date)"
samtools view -@ "$threads" -b -F 4 "$input_bam" > "$mapped_bam"

echo "[get_remap] intersecting mapped reads with centromeres: $(date)"
bedtools intersect -abam "$mapped_bam" -b "$centromere_bed" > "$centromere_bam"

echo "[get_remap] converting BAM subsets to FASTQ: $(date)"
samtools fastq -@ "$threads" "$unmapped_bam" > "$unmapped_fastq"
samtools fastq -@ "$threads" "$centromere_bam" > "$centromere_fastq"

echo "[get_remap] merging FASTQ subsets: $(date)"
cat "$centromere_fastq" "$unmapped_fastq" > "$merged_fastq"
gzip -c "$merged_fastq" > "$output_fastq_gz"

if [[ "$keep_intermediates" -eq 1 ]]; then
  cp "$unmapped_bam" "$outdir/${prefix}.unmapped.bam"
  cp "$mapped_bam" "$outdir/${prefix}.mapped.bam"
  cp "$centromere_bam" "$outdir/${prefix}.centromere.bam"
fi

echo "[get_remap] wrote $output_fastq_gz"
