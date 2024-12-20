#!/bin/bash

plink="/home/filippo/Downloads/plink"
basefolder="/home/filippo/Documents/ciampolini/golden_retriever"
inputfile="results/golden_filtered" ## binary ped files
outputfolder="results/ld"
max_snp_miss_rate=0.05
max_ind_miss_rate=0.10
min_maf=0.05
thin=0.90
ldwin=99999
ldwinr2=0.05
species=dog
option=autosome


if [[ ! -e ${basefolder}/${outputfolder} ]]; then
    mkdir ${basefolder}/${outputfolder}
fi

basefnm=${inputfile##*/}
echo "File name is $basefnm"

echo "##########################################"
echo " 1) Calculate LD:"
echo "##########################################"

## one chromosome for LD view
$plink --$species --bfile ${basefolder}/${inputfile} --r2 --maf $min_maf --ld-window $ldwin --ld-window-r2 $ldwinr2 --chr 38 --out ${basefolder}/${outputfolder}/${basefnm}_38

## all chromosomes
$plink --$species --bfile ${basefolder}/${inputfile} --r2 --maf $min_maf --ld-window $ldwin --ld-window-r2 $ldwinr2 --thin $thin --out ${basefolder}/${outputfolder}/${basefnm}


echo "DONE!!"
