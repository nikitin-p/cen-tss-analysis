#!/usr/bin/env bash

echo "STARTED at $(date)"

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
    "h9hesc1"
    "h9hesc2"
    "h9hesc3"
    "h9hesc4"
    "k56205M1"
    "k56205M2"
    "k5621M1"
    "k5621M2"
    "k5622M1"
    "k5622M2"
    "hepg205M1"
    "hepg205M2"
    "hepg21M1"
    "hepg21M2"
    "hepg22M1"
    "hepg22M2"
    "helas305M1"
    "helas305M2"
    "helas31M1"
    "helas31M2"
    "helas32M1"
    "helas32M2"
    "gm1287805M1"
    "gm1287805M2"
    "gm128781M1"
    "gm128781M2"
    "gm128782M1"
    "gm128782M2"
    "gm128782M3"
    "gm128782M4"
    "gm128784M1"
    "gm128784M2"
    "gm128786M1"
    "gm128786M2"
    "gm128788M1"
    "gm128788M2"
    "mcf705M1"
    "mcf705M2"
    "mcf71M1"
    "mcf71M2"
    "mcf72M1"
    "mcf72M2"
    "mcf72M3"
    "mcf72M4"
    "mcf74M1"
    "mcf74M2"
    "mcf76M1"
    "mcf76M2"
    "mcf78M1"
    "mcf78M2"
    "luhmesDay1Rep1NETCAGE"
    "luhmesDay1Rep2NETCAGE"
    "luhmesDay1Rep3NETCAGE"
    "luhmesDay3Rep1NETCAGE"
    "luhmesDay3Rep2NETCAGE"
    "luhmesDay3Rep3NETCAGE"
    "luhmesDay6Rep1NETCAGE"
    "luhmesDay6Rep2NETCAGE"
    "luhmesDay6Rep3NETCAGE"
)

for SAMPLE in "${SAMPLES[@]}"
do
    echo "Processing ${SAMPLE}"

    #samtools view -@ 20 -H ${SAMPLE}_bowtie2_centromeres.bam > ${SAMPLE}_header.BAM.sam

    #samtools view -@ 20 -F 20 ${SAMPLE}_bowtie2_centromeres.bam | awk -F '\t' 'BEGIN {OFS="\t"} {BASE = substr($10, 1, 1); if ($6 ~ /^1S[0-9]/ && BASE == "G") {print $0}}' > ${SAMPLE}_SoftclipG_F.BAM.sam

    #samtools view -@ 20 -f 16 ${SAMPLE}_bowtie2_centromeres.bam | awk -F '\t' 'BEGIN {OFS="\t"} { ALT = substr($10, length($10)-1, 1); if ($6 ~ /[0-9]M1S$/ && ALT == "C") {print $0}}' > ${SAMPLE}_SoftclipG_R.BAM.sam

    #cat ${SAMPLE}_header.BAM.sam ${SAMPLE}_SoftclipG_F.BAM.sam ${SAMPLE}_SoftclipG_R.BAM.sam | \
    #    samtools sort -@ 20 -O bam -o ${SAMPLE}.BAM_SoftclipG_centromeres.bam

    #rm ${SAMPLE}_header.BAM.sam
    #rm ${SAMPLE}_SoftclipG_F.BAM.sam
    #rm ${SAMPLE}_SoftclipG_R.BAM.sam

    samtools view -h ${SAMPLE}.BAM_SoftclipG_centromeres.bam | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {for (i=12; i<=NF; i++) if ($i ~ /^XM:i:/) {split($i, xm, ":"); if (xm[3]/length($10) < 0.10) print; break}}' | samtools view -b -o ${SAMPLE}.whole_centromere.1G_XM_90.bam

    samtools index ${SAMPLE}.whole_centromere.1G_XM_90.bam

    echo "DONE HISAT for ${SAMPLE} at $(date)"

done

echo "DONE at $(date)"