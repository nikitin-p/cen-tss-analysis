#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: bamcoverage_5p.sh -i INPUT_BAM -o OUTPUT_PREFIX [-t THREADS]

Generates two strand-specific 5'-end coverage bigWig files:
  OUTPUT_PREFIX.fwd.bw   forward strand (+)
  OUTPUT_PREFIX.rev.bw   reverse strand (-)

Both files use positive values. Coverage is counted at the 5'-end of each
alignment (--Offset 1 1 --binSize 1), which gives 1 bp resolution at the
first mapped position of each read.

IMPORTANT - 1 bp offset for initiator dinucleotide analysis:
  Reads in this pipeline carry a 5'-soft-clipped G (1S CIGAR prefix).
  bamCoverage reports the first MAPPED base, not the soft-clipped G.
  The actual G nucleotide is therefore:
    - forward strand: 1 bp upstream  of the reported bigWig position
    - reverse strand: 1 bp downstream of the reported bigWig position
  When extracting initiator dinucleotides in R, shift genomic coordinates
  by 1 accordingly (e.g. using rtracklayer + BSgenome).
USAGE
}

threads=8
input_bam=""
output_prefix=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-bam)      input_bam="$2";      shift 2 ;;
    -o|--output-prefix)  output_prefix="$2";  shift 2 ;;
    -t|--threads)        threads="$2";        shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$input_bam" ]]     || { echo "Missing --input-bam"     >&2; exit 1; }
[[ -n "$output_prefix" ]] || { echo "Missing --output-prefix" >&2; exit 1; }
[[ -f "$input_bam" ]]     || { echo "Input BAM not found: $input_bam" >&2; exit 1; }

mkdir -p "$(dirname "$output_prefix")"

echo "[bamcoverage_5p] forward strand: $(date)"
bamCoverage \
  -b "$input_bam" \
  -o "${output_prefix}.fwd.bw" \
  -p "$threads" \
  --Offset 1 1 \
  --binSize 1 \
  --samFlagExclude 16

echo "[bamcoverage_5p] reverse strand: $(date)"
bamCoverage \
  -b "$input_bam" \
  -o "${output_prefix}.rev.bw" \
  -p "$threads" \
  --Offset 1 1 \
  --binSize 1 \
  --samFlagInclude 16

echo "[bamcoverage_5p] wrote ${output_prefix}.fwd.bw"
echo "[bamcoverage_5p] wrote ${output_prefix}.rev.bw"
