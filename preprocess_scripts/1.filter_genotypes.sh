#!/bin/bash

plink="~/Downloads/plink"
basefolder="/home/diazjr/filippo/Montana"
inputfile="geno_montana_3808_51686" ## binary ped files
outputfolder="results"
max_snp_miss_rate=0.05

echo "filtering for SNP missing rate"
$plink --cow --bfile ${basefolder}/${inputfile} --geno $max_snp_miss_rate --make-bed --out "${basefolder}/${outputfolder}/montana_filtered"

echo "filtering for sample missing rate"
$plink --cow --bfile ${basefolder}/${outputfolder}/montana_filtered --geno $max_snp_miss_rate --make-bed --out "${basefolder}/${outputfolder}/montana_filtered"



# plink --cow --bfile geno_montana_3808_51686_geno05 --mind 0.05 --make-bed --out geno_montana_3808_51686_geno05_min05

# 1.4 Minor allele frequency ----------------------------------------------

# plink --cow --bfile geno_montana_3808_51686_geno05_min05 --freq –out geno_montana_3808_51686_geno05_min05_frq

echo "DONE!!"
