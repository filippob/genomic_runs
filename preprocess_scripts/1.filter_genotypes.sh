#!/bin/bash

plink="/home/filippo/Downloads/plink"
basefolder="/home/filippo/Documents/ciampolini/golden_retriever"
inputfile="dati/plink_files/golden" ## binary ped files
outputfolder="results"
max_snp_miss_rate=0.05
max_ind_miss_rate=0.10
min_maf=
species=dog
option=autosome-xy

if [[ ! -e ${basefolder}/${outputfolder} ]]; then
    mkdir ${basefolder}/${outputfolder}
fi

basefnm=${inputfile##*/}
echo "File name is $basefnm"

echo "##########################################"
echo " 1) Calculating missing rate before filtering:"
echo "##########################################"
$plink --$species --bfile ${basefolder}/${inputfile} --missing --out ${basefolder}/${outputfolder}/${basefnm}

echo "##########################################"
echo " 2) filtering for SNP missing rate"
echo "##########################################"
$plink --$species --bfile ${basefolder}/${inputfile} --geno $max_snp_miss_rate --$option --snps-only --make-bed --out "${basefolder}/${outputfolder}/${basefnm}_filtered"

echo "##########################################"
echo " 3) filtering for sample missing rate"
echo "##########################################"
$plink --$species --bfile ${basefolder}/${outputfolder}/${basefnm}_filtered --mind $max_ind_miss_rate --make-bed --out "${basefolder}/${outputfolder}/${basefnm}_filtered"

## MAF filtering: not needed for ROH analysis !!

if [ ! -z "${min_maf}" ]; then

	echo "##########################################"
	echo " 3.5) filtering for MAF"
	echo "##########################################"
	$plink --$species --bfile ${basefolder}/${outputfolder}/${basefnm}_filtered --maf $min_maf --make-bed --out ${basefolder}/${outputfolder}/${basefnm}_filtered
fi

echo "##########################################"
echo " 4) Calculating mssing rate after filtering:"
echo "##########################################"
$plink --$species --bfile ${basefolder}/${outputfolder}/${basefnm}_filtered --freq --missing --out ${basefolder}/${outputfolder}/${basefnm}_filtered

echo "DONE!!"
