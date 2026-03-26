#!/usr/bin/env bash

echo 'STARTED' at $(date)

SAMPLES=(
    "luhmesDay1Rep1CAGE"
    "luhmesDay1Rep2CAGE"
    "luhmesDay1Rep3CAGE"
    "luhmesDay3Rep1CAGE"
    "luhmesDay3Rep2CAGE"
    "luhmesDay3Rep3CAGE"
    "luhmesDay6Rep1CAGE"
    "luhmesDay6Rep2CAGE"
    "luhmesDay6Rep3CAGE"
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

    STAR \
        --runThreadN 20 \
        --genomeDir /home/pnikitin/cage/star \
        --readFilesIn /home/pnikitin/cage/data_net_cage/neuro/${SAMPLE}_L001_R1_001.fastq.gz \
        --outSAMtype BAM Unsorted \
        --readFilesCommand gunzip -c \
        --alignEndsType Local \
        --outSAMunmapped Within \
        --outFilterMultimapNmax 10 \
        --winAnchorMultimapNmax 20 \
        --outSAMmultNmax 1 \
        --outMultimapperOrder Random \
        --outFileNamePrefix ${SAMPLE}_initial_map_

    echo 'DONE STAR for' ${SAMPLE} $(date)
done