#!/bin/bash

plink="/home/filippo/Downloads/plink"
basefolder="/home/filippo/Documents/ciampolini/golden_retriever"
inputfile="dati/plink_files/golden" ## binary ped files
outputfolder="results"
max_snp_miss_rate=0.05
min_maf=0.05
species=dog
option=autosome-xy

if [[ ! -e ${basefolder}/${outputfolder} ]]; then
    mkdir ${basefolder}/${outputfolder}
fi

basefnm=${inputfile##*/}
echo "File name is $basefnm"

echo "filtering for SNP missing rate"
$plink --$species --bfile ${basefolder}/${inputfile} --geno $max_snp_miss_rate --$option --snps-only --make-bed --out "${basefolder}/${outputfolder}/${basefnm}_filtered"

echo "filtering for sample missing rate"
$plink --$species --bfile ${basefolder}/${outputfolder}/${basefnm}_filtered --mind $max_snp_miss_rate --make-bed --out "${basefolder}/${outputfolder}/${basefnm}_filtered"

## MAF filtering: not needed for ROH analysis !!
echo "filtering for MAF"
#$plink --$species --bfile ${basefolder}/${outputfolder}/${basefnm}_filtered --maf $min_maf --make-bed --out ${basefolder}/${outputfolder}/${basefnm}_filtered

# 1.4 Minor allele frequency ----------------------------------------------
$plink --$species --bfile ${basefolder}/${outputfolder}/${basefnm}_filtered --freq --missing --out ${basefolder}/${outputfolder}/${basefnm}_filtered

echo "DONE!!"
