#!/usr/bin/env bash

echo 'STARTED' at $(date)

for SAMPLE in ErythroidProgenitor HematopoieticProgenitor KeratinocyteProgenitor MyeloidProgenitor
do
    trim_galore \
        --fastqc \
        --clip_R1 9 \
        --three_prime_clip_R1 3 \
        --cores 10 \
        --gzip \
        ${SAMPLE}_L001_R1_001.fastq.gz
    echo 'DONE' ${SAMPLE} $(date)
done

for SAMPLE in FibroblastDermisAniso1 FibroblastDermisAniso2 OsteoblastAniso1 OsteoblastAniso2 SkeletalMuscleSatelliteAniso1 SkeletalMuscleSatelliteAniso2 VillousMesenchymeAniso1 VillousMesenchymeAniso2
do
    trim_galore \
        --fastqc \
        --clip_R1 9 \
        --three_prime_clip_R1 16 \
        --cores 10 \
        --gzip \
        ${SAMPLE}_L001_R1_001.fastq.gz
    echo 'DONE' ${SAMPLE} $(date)
done

echo 'STARTED' at $(date)

for SAMPLE in H9hESC
do
    trim_galore \
        --fastqc \
        --three_prime_clip_R1 2 \
        --cores 8 \
        --gzip \
        -o /home/pnikitin/cage/data/trimmed \
        ${SAMPLE}_L001_R1_001.fastq.gz
    echo 'DONE' ${SAMPLE} $(date)
done

for SAMPLE in KeratinocyteFemale UmbilicalVeinCytosolicFraction
do
    trim_galore \
        --fastqc \
        --clip_R1 9 \
        --three_prime_clip_R1 1 \
        --cores 8 \
        --gzip \
        -o /home/pnikitin/cage/data/trimmed \
        ${SAMPLE}_L001_R1_001.fastq.gz
    echo 'DONE' ${SAMPLE} $(date)
done

for SAMPLE in BCell CD14PositiveMonocyte ArticularChondrocyteKneeJoint BronchialSmoothMuscleCells
do
    trim_galore \
        --fastqc \
        --clip_R1 9 \
        --three_prime_clip_R1 16 \
        --cores 8 \
        --gzip \
        -o /home/pnikitin/cage/data/trimmed \
        ${SAMPLE}_L001_R1_001.fastq.gz
    echo 'DONE' ${SAMPLE} $(date)
done

echo 'FINISHED' at $(date)