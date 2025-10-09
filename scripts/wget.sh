#!/usr/bin/env bash

# Keratinocyte progenitor SE
# MD5 046ff503b4ed35daba076fdd91fea607
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/003/SRR1723253/SRR1723253.fastq.gz \
    -O keratinocyte_progenitor.fastq.gz

# Hematopoietic stem/progenitor cells SE
# MD5 046d04d155d0b1e6e3c78019420fd88f
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/000/SRR1720190/SRR1720190.fastq.gz \
    -O hematopoietic_progenitor.fastq.gz

# Erythroid progenitors/precursors SE
# MD5 095915bd53040a3f59e556f7ebd3c75c
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/001/SRR1720191/SRR1720191.fastq.gz \
    -O erythroid_progenitor.fastq.gz

# Myeloid progenitors/precursors SE
# MD5 d728a0340c8145e19dc60e368690b3ee
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/002/SRR1720192/SRR1720192.fastq.gz \
    -O myeloid_progenitor.fastq.gz

# Fibroblast of dermis aniso1 SE
# MD5 b43d322b489b0d416e7c2014cd356a47
wget \
    -nc https://www.encodeproject.org/files/ENCFF000USW/@@download/ENCFF000USW.fastq.gz \
    -O fibroblast_dermis_aniso1.fastq.gz

# Fibroblast of dermis aniso2 SE
# MD5 7c15e090180a3b883344bf0f2ac3d730
wget \
    -nc https://www.encodeproject.org/files/ENCFF000USV/@@download/ENCFF000USV.fastq.gz \
    -O fibroblast_dermis_aniso2.fastq.gz

# Skeletal muscle satellite cell aniso1 SE
# MD5 2e86d7546804759d5105ce1d9deb9bc6
wget \
    -nc https://www.encodeproject.org/files/ENCFF000UWW/@@download/ENCFF000UWW.fastq.gz \
    -O skeletal_muscle_satellite_aniso1.fastq.gz

# Skeletal muscle satellite cell aniso2 SE
# MD5 7bc5587bcdd437ed5562fe89d9ece389
wget \
    -nc https://www.encodeproject.org/files/ENCFF000UWZ/@@download/ENCFF000UWZ.fastq.gz \
    -O skeletal_muscle_satellite_aniso2.fastq.gz

# Keratinocyte iso1 SE
# MD5 e36c75543fd9e57d3254ab9ab082d6be
wget \
    -nc https://www.encodeproject.org/files/ENCFF000UTO/@@download/ENCFF000UTO.fastq.gz \
    -O keratinocyte_iso1.fastq.gz

# Keratinocyte iso2 SE
# MD5 4e735ecd76e55d3e1b47f0b299b67adf
wget \
    -nc https://www.encodeproject.org/files/ENCFF000UTQ/@@download/ENCFF000UTQ.fastq.gz \
    -O keratinocyte_iso2.fastq.gz

# Villous Mesenchyme Aniso1 SE 51
# MD5 5fc59df23edea28a3d417e10d474fe45
wget \
    -nc https://www.encodeproject.org/files/ENCFF000UKV/@@download/ENCFF000UKV.fastq.gz \
    -O VillousMesenchymeAniso1_L001_R1_001.fastq.gz

# Villous Mesenchyme Aniso2 SE 51
# MD5 101789c4fb73ae904a468529be5063f0
wget \
    -nc https://www.encodeproject.org/files/ENCFF000ULB/@@download/ENCFF000ULB.fastq.gz \
    -O VillousMesenchymeAniso2_L001_R1_001.fastq.gz

# Osteoblast Aniso1 SE 51
# MD5 b2b468f53cfffef591e5373590e2baf2
wget \
    -nc https://www.encodeproject.org/files/ENCFF000UHD/@@download/ENCFF000UHD.fastq.gz \
    -O OsteoblastAniso1_L001_R1_001.fastq.gz

# Osteoblast Aniso2 SE 51
# MD5 f0eb0e43c546394631334cfcab12379d
wget \
    -nc https://www.encodeproject.org/files/ENCFF000UHH/@@download/ENCFF000UHH.fastq.gz \
    -O OsteoblastAniso2_L001_R1_001.fastq.gz

# h9hesc1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/027/SRR14199827/SRR14199827.fastq.gz \
    -O h9hesc1_L001_R1_001.fastq.gz
# h9hesc2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/028/SRR14199828/SRR14199828.fastq.gz \
    -O h9hesc2_L001_R1_001.fastq.gz
# h9hesc3
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/029/SRR14199829/SRR14199829.fastq.gz \
    -O h9hesc3_L001_R1_001.fastq.gz
# h9hesc4
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/030/SRR14199830/SRR14199830.fastq.gz \
    -O h9hesc4_L001_R1_001.fastq.gz
# k56205M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/000/SRR7641290/SRR7641290.fastq.gz \
    -O k56205M1_L001_R1_001.fastq.gz
# k56205M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/003/SRR7641293/SRR7641293.fastq.gz \
    -O k56205M2_L001_R1_001.fastq.gz
# k5621M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/001/SRR7641291/SRR7641291.fastq.gz \
    -O k5621M1_L001_R1_001.fastq.gz
# k5621M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/004/SRR7641294/SRR7641294.fastq.gz \
    -O k5621M2_L001_R1_001.fastq.gz
# k5622M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/002/SRR7641292/SRR7641292.fastq.gz \
    -O k5622M1_L001_R1_001.fastq.gz
# k5622M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/005/SRR7641295/SRR7641295.fastq.gz \
    -O k5622M2_L001_R1_001.fastq.gz
# hepg205M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/002/SRR7641282/SRR7641282.fastq.gz \
    -O hepg205M1_L001_R1_001.fastq.gz
# hepg205M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/005/SRR7641285/SRR7641285.fastq.gz \
    -O hepg205M2_L001_R1_001.fastq.gz
# hepg21M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/003/SRR7641283/SRR7641283.fastq.gz \
    -O hepg21M1_L001_R1_001.fastq.gz
# hepg21M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/006/SRR7641286/SRR7641286.fastq.gz \
    -O hepg21M2_L001_R1_001.fastq.gz
# hepg22M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/004/SRR7641284/SRR7641284.fastq.gz \
    -O hepg22M1_L001_R1_001.fastq.gz
# hepg22M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/007/SRR7641287/SRR7641287.fastq.gz \
    -O hepg22M2_L001_R1_001.fastq.gz
# helas305M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/004/SRR7641274/SRR7641274.fastq.gz \
    -O helas305M1_L001_R1_001.fastq.gz
# helas305M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/007/SRR7641277/SRR7641277.fastq.gz \
    -O helas305M2_L001_R1_001.fastq.gz
# helas31M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/005/SRR7641275/SRR7641275.fastq.gz \
    -O helas31M1_L001_R1_001.fastq.gz
# helas31M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/008/SRR7641278/SRR7641278.fastq.gz \
    -O helas31M2_L001_R1_001.fastq.gz
# helas32M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/006/SRR7641276/SRR7641276.fastq.gz \
    -O helas32M1_L001_R1_001.fastq.gz
# helas32M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/009/SRR7641279/SRR7641279.fastq.gz \
    -O helas32M2_L001_R1_001.fastq.gz
# gm1287805M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/008/SRR7641258/SRR7641258.fastq.gz \
    -O gm1287805M1_L001_R1_001.fastq.gz
# gm1287805M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/001/SRR7641261/SRR7641261.fastq.gz \
    -O gm1287805M2_L001_R1_001.fastq.gz
# gm128781M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/009/SRR7641259/SRR7641259.fastq.gz \
    -O gm128781M1_L001_R1_001.fastq.gz
# gm128781M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/002/SRR7641262/SRR7641262.fastq.gz \
    -O gm128781M2_L001_R1_001.fastq.gz
# gm128782M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/000/SRR7641260/SRR7641260.fastq.gz \
    -O gm128782M1_L001_R1_001.fastq.gz
# gm128782M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/003/SRR7641263/SRR7641263.fastq.gz \
    -O gm128782M2_L001_R1_001.fastq.gz
# gm128782M3
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/004/SRR7641264/SRR7641264.fastq.gz \
    -O gm128782M3_L001_R1_001.fastq.gz
# gm128782M4
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/008/SRR7641268/SRR7641268.fastq.gz \
    -O gm128782M4_L001_R1_001.fastq.gz
# mcf705M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/005/SRR7641185/SRR7641185.fastq.gz \
    -O mcf705M1_L001_R1_001.fastq.gz
# mcf705M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/008/SRR7641188/SRR7641188.fastq.gz \
    -O mcf705M2_L001_R1_001.fastq.gz
# mcf71M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/006/SRR7641186/SRR7641186.fastq.gz \
    -O mcf71M1_L001_R1_001.fastq.gz
# mcf71M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/009/SRR7641189/SRR7641189.fastq.gz \
    -O mcf71M2_L001_R1_001.fastq.gz
# mcf72M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/007/SRR7641187/SRR7641187.fastq.gz \
    -O mcf72M1_L001_R1_001.fastq.gz
# mcf72M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/000/SRR7641190/SRR7641190.fastq.gz \
    -O mcf72M2_L001_R1_001.fastq.gz
# mcf72M3
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/001/SRR7641191/SRR7641191.fastq.gz \
    -O mcf72M3_L001_R1_001.fastq.gz
# mcf72M4
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/005/SRR7641195/SRR7641195.fastq.gz \
    -O mcf72M4_L001_R1_001.fastq.gz

# h9hescCAGE1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/019/SRR14199819/SRR14199819.fastq.gz \
    -O h9hescCAGE1_L001_R1_001.fastq.gz
# h9hescCAGE2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/020/SRR14199820/SRR14199820.fastq.gz \
    -O h9hescCAGE2_L001_R1_001.fastq.gz
# h9hescCAGE3
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/021/SRR14199821/SRR14199821.fastq.gz \
    -O h9hescCAGE3_L001_R1_001.fastq.gz
# h9hescCAGE4
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/022/SRR14199822/SRR14199822.fastq.gz \
    -O h9hescCAGE4_L001_R1_001.fastq.gz
# k562CAGE1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/008/SRR7641288/SRR7641288.fastq.gz \
    -O k562CAGE1_L001_R1_001.fastq.gz
# k562CAGE2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/009/SRR7641289/SRR7641289.fastq.gz \
    -O k562CAGE2_L001_R1_001.fastq.gz
# hepg2CAGE1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/000/SRR7641280/SRR7641280.fastq.gz \
    -O hepg2CAGE1_L001_R1_001.fastq.gz
# hepg2CAGE2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/001/SRR7641281/SRR7641281.fastq.gz \
    -O hepg2CAGE2_L001_R1_001.fastq.gz
# helas3CAGE1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/002/SRR7641272/SRR7641272.fastq.gz \
    -O helas3CAGE1_L001_R1_001.fastq.gz
# helas3CAGE2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/003/SRR7641273/SRR7641273.fastq.gz \
    -O helas3CAGE2_L001_R1_001.fastq.gz
# gm12878CAGE1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/006/SRR7641256/SRR7641256.fastq.gz \
    -O gm12878CAGE1_L001_R1_001.fastq.gz
# gm12878CAGE2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/007/SRR7641257/SRR7641257.fastq.gz \
    -O gm12878CAGE2_L001_R1_001.fastq.gz
# mcf7CAGE1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/003/SRR7641183/SRR7641183.fastq.gz \
    -O mcf7CAGE1_L001_R1_001.fastq.gz
# mcf7CAGE2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/004/SRR7641184/SRR7641184.fastq.gz \
    -O mcf7CAGE2_L001_R1_001.fastq.gz

# gm128784M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/005/SRR7641265/SRR7641265.fastq.gz \
    -O gm128784M1_L001_R1_001.fastq.gz
# gm128784M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/009/SRR7641269/SRR7641269.fastq.gz \
    -O gm128784M2_L001_R1_001.fastq.gz
# gm128786M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/006/SRR7641266/SRR7641266.fastq.gz \
    -O gm128786M1_L001_R1_001.fastq.gz
# gm128786M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/000/SRR7641270/SRR7641270.fastq.gz \
    -O gm128786M2_L001_R1_001.fastq.gz
# gm128788M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/007/SRR7641267/SRR7641267.fastq.gz \
    -O gm128788M1_L001_R1_001.fastq.gz
# gm128788M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/001/SRR7641271/SRR7641271.fastq.gz \
    -O gm128788M2_L001_R1_001.fastq.gz
# mcf74M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/002/SRR7641192/SRR7641192.fastq.gz \
    -O mcf74M1_L001_R1_001.fastq.gz
# mcf74M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/006/SRR7641196/SRR7641196.fastq.gz \
    -O mcf74M2_L001_R1_001.fastq.gz
# mcf746M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/003/SRR7641193/SRR7641193.fastq.gz \
    -O mcf746M1_L001_R1_001.fastq.gz
# mcf76M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/007/SRR7641197/SRR7641197.fastq.gz \
    -O mcf76M2_L001_R1_001.fastq.gz
# mcf78M1
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/004/SRR7641194/SRR7641194.fastq.gz \
    -O mcf78M1_L001_R1_001.fastq.gz
# mcf78M2
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/008/SRR7641198/SRR7641198.fastq.gz \
    -O mcf78M2_L001_R1_001.fastq.gz



# luhmesDay1Rep1CAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/098/ERR12322798/ERR12322798.fastq.gz \
    -O luhmesDay1Rep1CAGE_L001_R1_001.fastq.gz
# luhmesDay1Rep2CAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/099/ERR12322799/ERR12322799.fastq.gz \
    -O luhmesDay1Rep2CAGE_L001_R1_001.fastq.gz
# luhmesDay1Rep3CAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/000/ERR12322800/ERR12322800.fastq.gz \
    -O luhmesDay1Rep3CAGE_L001_R1_001.fastq.gz
# luhmesDay3Rep1CAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/001/ERR12322801/ERR12322801.fastq.gz \
    -O luhmesDay3Rep1CAGE_L001_R1_001.fastq.gz
# luhmesDay3Rep2CAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/002/ERR12322802/ERR12322802.fastq.gz \
    -O luhmesDay3Rep2CAGE_L001_R1_001.fastq.gz
# luhmesDay3Rep3CAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/003/ERR12322803/ERR12322803.fastq.gz \
    -O luhmesDay3Rep3CAGE_L001_R1_001.fastq.gz
# luhmesDay6Rep1CAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/004/ERR12322804/ERR12322804.fastq.gz \
    -O luhmesDay6Rep1CAGE_L001_R1_001.fastq.gz
# luhmesDay6Rep2CAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/005/ERR12322805/ERR12322805.fastq.gz \
    -O luhmesDay6Rep2CAGE_L001_R1_001.fastq.gz
# another try
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR123/ERR12322805/Total.Day6.S06.fastq.bz2 \
    -O luhmesDay6Rep2CAGE_L001_R1_001.fastq.gz
# luhmesDay6Rep3CAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/006/ERR12322806/ERR12322806.fastq.gz \
    -O luhmesDay6Rep3CAGE_L001_R1_001.fastq.gz


# luhmesDay1Rep1NETCAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/089/ERR12322789/ERR12322789.fastq.gz \
    -O luhmesDay1Rep1NETCAGE_L001_R1_001.fastq.gz
# luhmesDay1Rep2NETCAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/090/ERR12322790/ERR12322790.fastq.gz \
    -O luhmesDay1Rep2NETCAGE_L001_R1_001.fastq.gz
# luhmesDay1Rep3NETCAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/091/ERR12322791/ERR12322791.fastq.gz \
    -O luhmesDay1Rep3NETCAGE_L001_R1_001.fastq.gz
# luhmesDay3Rep1NETCAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/092/ERR12322792/ERR12322792.fastq.gz \
    -O luhmesDay3Rep1NETCAGE_L001_R1_001.fastq.gz
# luhmesDay3Rep2NETCAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/093/ERR12322793/ERR12322793.fastq.gz \
    -O luhmesDay3Rep2NETCAGE_L001_R1_001.fastq.gz
# luhmesDay3Rep3NETCAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/094/ERR12322794/ERR12322794.fastq.gz \
    -O luhmesDay3Rep3NETCAGE_L001_R1_001.fastq.gz
# luhmesDay6Rep1NETCAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/095/ERR12322795/ERR12322795.fastq.gz \
    -O luhmesDay6Rep1NETCAGE_L001_R1_001.fastq.gz
# luhmesDay6Rep2NETCAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/096/ERR12322796/ERR12322796.fastq.gz \
    -O luhmesDay6Rep2NETCAGE_L001_R1_001.fastq.gz
# luhmesDay6Rep3NETCAGE
wget \
    -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/097/ERR12322797/ERR12322797.fastq.gz \
    -O luhmesDay6Rep3NETCAGE_L001_R1_001.fastq.gz
