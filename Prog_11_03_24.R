#11_03_2024
# Programacion - 
# Poblacion Colombiana de Ovejas

#The input files for detectRUNS are Plink ped/map files

# Guidelines Genomic characterization of animal genetic resources FAO
#https://www.fao.org/documents/card/en/c/cc3079en 
#missingness
## plink --sheep --allow-no-sex --allow-extra-chr --bfile Colombia_updated_name_chr_bp 
##    --no-pheno --autosome --make-bed --missing --out missing_sheep_colombia

library(dplyr)
library(data.table)
# 1. Quality_Control ---------------------------------------------------------
# ## 1.1. Missingness -----------------------------------------------------------
library(data.table)

imiss<- fread("missing_sheep_colombia.imiss")
lmiss <- fread("missing_sheep_colombia.lmiss")
par(mar = c(2.5, 4.5, 1.5, 1.5))
par(mfrow =c(1,1))
hist(imiss$F_MISS, 
     xlab="missigness freq", cex=0.5,
     main="Individual missigness",cex.axis=0.7, las = 1,  
     breaks=50)

par(mar = c(2.5, 4.5, 1.5, 1.5))
hist(lmiss$F_MISS, 
     xlab="missigness freq", cex=0.5,
     main="Loci missigness",cex.axis=0.7, las = 1,    # Horizontal labels
     getOption("scipen"),
     breaks=50)

# To maximize the number of individuals, first prune for loci missingness, followed by individual missingness. Pruning for missingness can be performed in PLINK using the --geno
#and --mind for loci and individuals, respectively.
######  geno 0.05
#  plink --sheep --bfile Colombia_updated_name_chr_bp --geno 0.05  --make-bed -- out colombia_geno05



#mind 0.05
#plink --sheep --bfile colombia_geno05 --mind 0.05 --make-bed --out colombia_geno05_mind05

# The missingness statistics are then checked again:
# plink --sheep --bfile colombia_geno05_mind05 --missing --out colombia_geno05_mind05

imiss_after<- fread("colombia_geno05_mind05.imiss")
lmiss_after <- fread("colombia_geno05_mind05.lmiss")
par(mar = c(2.5, 4.5, 1.5, 1.5))
par(mfrow =c(1,1))
hist(imiss_after$F_MISS, 
     xlab="missigness freq", cex=0.5,
     main="Individual missigness-after pruning",cex.axis=0.7, las = 1,  
     col ="green",
     breaks=50)
par(mar = c(2.5, 4.5, 1.5, 1.5))
hist(lmiss_after$F_MISS, 
     xlab="missigness freq", cex=0.5,
     main="Loci missigness -after pruning",cex.axis=0.7, las = 1,    # Horizontal labels
     col ="green",
     getOption("scipen"),
     breaks=50)

#plink --sheep --bfile colombia_geno05_mind05 --maf 0.01 --make-bed --out colombia_geno05_mind05_maf001

genotypeFilePath <- "colombia_geno05_mind05_maf001.ped"
mapFilePath <- "colombia_geno05_mind05_maf001.map"

# genotypeFilePath <- "Brasil.ped"
# mapFilePath <- "Brasil.map"

slidingRuns_het <- slidingRUNS.run(
  genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath, 
  windowSize = 15, 
  threshold = 0.05,
  minSNP = 10, 
  ROHet = TRUE, 
  maxOppWindow = 1, 
  maxMissWindow = 1,
  maxGap = 10^6, 
  minLengthBps = 10000, 
  minDensity = 0.1, # SNP/kbps
  maxOppRun = 2,
  maxMissRun = 1
)

consecutiveRuns_het <- consecutiveRUNS.run(
  genotypeFile =genotypeFilePath,
  mapFile = mapFilePath,
  minSNP = 10,
  ROHet = TRUE,
  maxGap = 10^6,
  minLengthBps = 10000,
  maxOppRun = 2,
  maxMissRun = 1
)

summaryList <- summaryRuns(
  runs = consecutiveRuns_het, mapFile = mapFilePath, genotypeFile = genotypeFilePath, 
  Class = 6, snpInRuns = TRUE)
summaryList$summary_ROH_count
summaryList$summary_ROH_mean_chr
head(summaryList$SNPinRun)
plot_Runs(runs = consecutiveRuns_het)
plot_StackedRuns(runs = consecutiveRuns_het)
plot_SnpsInRuns(
  runs = consecutiveRuns_het[consecutiveRuns_het$chrom==2,], genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath)

plot_manhattanRuns(
  runs = consecutiveRuns_het[consecutiveRuns_het$group=="Jacobs",], 
  genotypeFile = genotypeFilePath, mapFile = mapFilePath)
