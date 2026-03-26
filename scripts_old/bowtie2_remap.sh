#!/usr/bin/env bash

echo 'STARTED' at $(date)

SAMPLES=(
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

    bowtie2 \
        -k 800 \
        -D 1600 \
        -R 1 \
        --local \
        -x /home/pnikitin/cage/bowtie2/t2t_hor_with_strand_centromeres_ref \
        -U /home/pnikitin/cage/bam_star_bowtie/${SAMPLE}.BAM1_remap_R1.fastq.gz \
        -p 20 | \
        samtools view -b - | \
        samtools sort -@ 20 -O bam -o ${SAMPLE}_bowtie2_centromeres.bam

    echo 'DONE BOWTIE2 for' ${SAMPLE} $(date)
done