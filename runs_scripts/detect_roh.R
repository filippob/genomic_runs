# library("devtools")
# devtools::install_github("bioinformatics-ptp/detectRUNS/detectRUNS")
library("tidyverse")
library("detectRUNS")
library("data.table")

args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1) {
  
  #loading the parameters
  source(args[1])
  # source("Analysis/hrr/config.R")
  
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    base_folder = '/home/filippo/Documents/ciampolini/golden_retriever',
    genotypes = "results/roh/golden_filtered_thin.ped",
    mapfile = "results/roh/golden_filtered_thin.map",
    windowSize = 50,
    threshold = 0.05,
    minSNP = 20,
    ROHet = FALSE,
    maxGap = 2.5e+4,
    minLengthBps = 1e+6,
    maxOppRun = 0,
    maxMissRun = 0,
    minDensity = 1/(5e+4),
    maxMissWindow = 0,
    maxOppWindow = 2,
    sliding = FALSE,
    force_overwrite = FALSE
  ))
}

writeLines(' - input files')
print(paste("genotype file:", config$genotypes))
print(paste("mapfile:", config$mapfile))

## detect ROH

if (config$sliding) {
  
  writeLines(' - detecting ROH with the sliding method')
  
  writeLines(' - current values of parameters')
  print(paste("min n. of SNP:", config$minSNP))
  print(paste("maximum gap between SNPs:", config$maxGap))
  print(paste("min length of roh (bps):", config$minLengthBps))
  print(paste("max homozygous SNP on roh:", config$maxOppRun))
  print(paste("max missing SNP in roh:", config$maxMissRun))
  print(paste("min SNP density in sliding window:", config$minDensity))
  print(paste("max missing SNP in sliding window:", config$maxMissWindow))
  print(paste("max heterozygous SNP in sliding window:", config$maxOppWindow))
  print(paste("size of sliding window:", config$windowSize))
  print(paste("threshold to call sliding window in roh:", config$threshold))
  
  resultRuns <- slidingRUNS.run(
    genotypeFile = config$genotypes,
    mapFile = config$mapfile,
    windowSize = config$windowSize,
    threshold = config$threshold,
    minSNP = config$minSNP,
    ROHet = config$ROHet,
    maxGap = config$maxGap,
    minLengthBps = config$minLengthBps,
    maxOppRun = config$maxOppRun,
    maxMissRun = config$maxMissRun,
    minDensity = config$minDensity,
    maxMissWindow = config$maxMissWindow,
    maxOppWindow = config$maxOppWindow
  )
  
} else {
  
  writeLines(' - detecting ROH with the consecutive method (Marras et al., 2016)')
  
  writeLines(' - current values of parameters')
  print(paste("min n. of SNP:", config$minSNP))
  print(paste("maximum gap between SNPs:", config$maxGap))
  print(paste("min length of roh (bps):", config$minLengthBps))
  print(paste("max homozygous SNP on roh:", config$maxOppRun))
  print(paste("max missing SNP in roh:", config$maxMissRun))
  
  resultRuns <- consecutiveRUNS.run(
    genotypeFile = config$genotypes,
    mapFile = config$mapfile,
    minSNP = config$minSNP,
    ROHet = config$ROHet,
    maxGap = config$maxGap,
    minLengthBps = config$minLengthBps,
    maxOppRun = config$maxOppRun,
    maxMissRun = config$maxMissRun
  )
  
}

fname = file.path(config$base_folder, "results/roh/roh.csv")
writeLines(" - writing out results to file")
print(paste("output file name:",fname))
fwrite(resultRuns, file = fname, sep = ",")

