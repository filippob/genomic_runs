#!/bin/bash

## script to detect ROH/HRR
## the bash scripts first thins and transform the input data
## then calls an R script to detect HRR (detectRUNS inside)

plink="/home/filippo/Downloads/plink"
basefolder="/home/filippo/Documents/ciampolini/golden_retriever"
inputfile="results/golden_filtered" ## binary ped files
outputfolder="results/hrr"
species=dog
thin=0.99
minMAF=0.05
rscript=$HOME/Documents/SMARTER/hrr_goats/scripts/detect_hrr.R
minSNP=15
maxGap=10^6
minBps=250*10^3
maxOpp=3
maxMiss=0


if [ ! -d "${basefolder}/${outputfolder}" ]; then
	mkdir ${basefolder}/${outputfolder}
fi

basefnm=${inputfile##*/}
echo "File name is $basefnm"

## Thinning
echo "1. Thinning"

## plink commands 
## (option --cow to allow for 60 chromosomes in goats)
## there is no option --goat 
## (sheep has 54 chromosomes, therefore --sheep would not work)
$plink --$species --bfile ${basefolder}/${inputfile} --maf $minMAF --thin $thin --make-bed --out ${basefolder}/${outputfolder}/${basefnm}_thin

## when FID == IID, this may give problems with detectRUNS (especially plot functions)
echo "2. Update population information"

cut -f1-2 -d' ' ${basefolder}/${outputfolder}/${basefnm}_thin.fam > ${basefolder}/${outputfolder}/temp1
sed -i 's/\s*$/ GOLD/' ${basefolder}/${outputfolder}/temp1
cut -f1 -d' ' ${basefolder}/${outputfolder}/temp1 > ${basefolder}/${outputfolder}/temp2
paste ${basefolder}/${outputfolder}/temp1 ${basefolder}/${outputfolder}/temp2 > ${basefolder}/${outputfolder}/newpop.txt

$plink --$species --bfile ${basefolder}/${outputfolder}/${basefnm}_thin --update-ids ${basefolder}/${outputfolder}/newpop.txt --recode --out ${basefolder}/${outputfolder}/${basefnm}_thin

## Configuration file
echo "2. Creating the configuration file"

OUTDIR=${basefolder}/${outputfolder}

echo "config = data.frame(" > $OUTDIR/config.R
echo "base_folder = '~/Documents/SMARTER/Analysis/hrr/'," >> $OUTDIR/config.R
echo "genotypes = '$OUTDIR/goat_thin.ped'," >> $OUTDIR/config.R
echo "mapfile = '$OUTDIR/goat_thin.map'," >> $OUTDIR/config.R
echo "minSNP = $minSNP," >> $OUTDIR/config.R
echo "ROHet = TRUE," >> $OUTDIR/config.R
echo "maxGap = $maxGap," >> $OUTDIR/config.R
echo "minLengthBps = $minBps," >> $OUTDIR/config.R
echo "maxOppRun = $maxOpp," >> $OUTDIR/config.R
echo "maxMissRun = $maxMiss," >> $OUTDIR/config.R
echo "force_overwrite = FALSE)" >> $OUTDIR/config.R

## detect RUNS
echo "3. Detecting HRR"
#Rscript --vanilla $RSCRIPT $OUTDIR/config.R

## house cleaning
echo "4. Cleaning"
#rm $OUTDIR/goat_thin.log
#rm $OUTDIR/goat_thin.nosex
#rm $OUTDIR/config.R

echo "DONE!"
