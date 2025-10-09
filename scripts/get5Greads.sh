#!/usr/bin/env bash

echo 'STARTED' at $(date)

for SAMPLE in K562 HepG2 HeLaS3 GM12878 MCF7 FANTOM6DermalLymphaticVascularEndothelialCells FANTOM6DermalBloodVascularEndothelialCells FANTOM6DermalFibroblasts
#for SAMPLE in ErythroidProgenitor FibroblastDermisAniso1 FibroblastDermisAniso2 HematopoieticProgenitor KeratinocyteProgenitor MyeloidProgenitor OsteoblastAniso1 OsteoblastAniso2 SkeletalMuscleSatelliteAniso1 SkeletalMuscleSatelliteAniso2 VillousMesenchymeAniso1 VillousMesenchymeAniso2
#for SAMPLE in KeratinocyteFemale UmbilicalVeinCytosolicFraction K562 HepG2 HeLaS3 GM12878 MCF7 FANTOM6DermalLymphaticVascularEndothelialCells FANTOM6DermalBloodVascularEndothelialCells FANTOM6DermalFibroblasts ErythroidProgenitor FibroblastDermisAniso2 HematopoieticProgenitor KeratinocyteProgenitor MyeloidProgenitor OsteoblastAniso1 OsteoblastAniso2 SkeletalMuscleSatelliteAniso1 VillousMesenchymeAniso1
do
    zcat ${SAMPLE}_L001_R1_001.fastq.gz | awk '(NR % 4 == 2 && substr($0, 1, 1) == "G"){print last; print; getline; print; getline; print} {last=$0}' > ${SAMPLE}_L001_R1_001_G.fastq
    gzip ${SAMPLE}_L001_R1_001_G.fastq
    echo 'DONE' ${SAMPLE} $(date)
done