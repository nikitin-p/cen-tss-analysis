#!/usr/bin/env bash

echo 'STARTED' at $(date)

SAMPLES=(
    "h9hescCAGE1"
    "h9hescCAGE2"
    "h9hescCAGE3"
    "h9hescCAGE4"
    "k562CAGE1"
    "k562CAGE2"
    "hepg2CAGE1"
    "hepg2CAGE2"
    "helas3CAGE1"
    "helas3CAGE2"
    "gm12878CAGE1"
    "gm12878CAGE2"
    "mcf7CAGE1"
    "mcf7CAGE2"
    "luhmesDay1Rep1CAGE"
    "luhmesDay1Rep2CAGE"
    "luhmesDay1Rep3CAGE"
    "luhmesDay3Rep1CAGE"
    "luhmesDay3Rep2CAGE"
    "luhmesDay3Rep3CAGE"
    "luhmesDay6Rep1CAGE"
    "luhmesDay6Rep2CAGE"
    "luhmesDay6Rep3CAGE"
)

for SAMPLE in "${SAMPLES[@]}"
do
    echo ${SAMPLE}
    samtools view -@ 20 -H ${SAMPLE}_initial_map_Aligned.out.bam > ${SAMPLE}_header.BAM.sam

    samtools view -@ 20 -F 20 ${SAMPLE}_initial_map_Aligned.out.bam | \
        awk -F '\t' 'BEGIN {OFS="\t"} {
            if ($6 ~ /^1S[0-9]+M/) { print $0 }}' > ${SAMPLE}_SoftclipG_F.BAM.sam

    samtools view -@ 20 -f 16 ${SAMPLE}_initial_map_Aligned.out.bam | \
        awk -F '\t' 'BEGIN {OFS="\t"} {
            if ($6 ~ /[0-9]M1S$/) { print $0 }}' > ${SAMPLE}_SoftclipG_R.BAM.sam

    cat ${SAMPLE}_header.BAM.sam ${SAMPLE}_SoftclipG_F.BAM.sam ${SAMPLE}_SoftclipG_R.BAM.sam | \
        samtools sort -@ 20 -O bam -o ${SAMPLE}.BAM_SoftclipG_initial_map.bam

    rm ${SAMPLE}_header.BAM.sam
    rm ${SAMPLE}_SoftclipG_F.BAM.sam
    rm ${SAMPLE}_SoftclipG_R.BAM.sam

    echo "DONE for ${SAMPLE} at $(date)"
done