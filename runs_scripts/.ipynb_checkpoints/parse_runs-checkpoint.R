devtools::install_github("bioinformatics-ptp/detectRUNS/detectRUNS")

library("tidyverse")
library("detectRUNS")
library("data.table")
# library("devtools")
# devtools::install_github("bioinformatics-ptp/detectRUNS/detectRUNS", ref = "master")


### CONFIG FILE
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
    repo_name = "genomic_runs", ## name of git repository from github /filippob/heterozygosity_rich_regions
    runs_file = "results/roh/roh.csv",
    pedfile = "results/roh/golden_filtered_thin.ped",
    mapfile = "results/roh/golden_filtered_thin.map",
    outdir = "results/roh",
    prefix = "roh_sliding", ## analysis identifier,
    thr_value = 0.20, ##threshold for ROH/HRR islands
    force_overwrite = FALSE
  ))
}

writeLines(' - current values of parameters')
print(paste("base folder:", config$base_folder))
print(paste("runs file:", config$runs_file))
print(paste("pedfile:", config$pedfile))
print(paste("mapfile:", config$mapfile))
print(paste("output folder:", config$outdir))
print(paste("repository name:", config$repo_name))
print(paste("identifier:", config$prefix))
print(paste("HRR threshold:", config$thr_value))

### READ DATA
writeLines(' - reading data')
runs = fread(config$runs_file)
pedfilePath = config$pedfile
mapfilePath = config$mapfile

## ROH/HRR summary stats
writeLines(' - calculating ROH/HRR descriptive stats')
fname = paste(config$prefix, "_summary_runs.csv", sep="")
fname = file.path(config$base_folder, config$outdir, fname)

runs %>%
  group_by(group) %>%
  summarise(N=n(), n_ind = length(unique(id)), avg_n = N/n_ind, avg_length = mean(lengthBps)) %>%
  fwrite(fname, sep = ",", col.names = TRUE)

## plot runs per samples
writeLines(' - plotting individual runs')
# runs$group <- factor(runs$group, levels = c("ALP","BOE","SAA","BRK","CRE","LNR")) ## per-breed data
fname = paste(config$prefix, "_plot_runs.pdf", sep="")
fname = file.path(config$base_folder, config$outdir, fname)
plot_Runs(runs = mutate(runs, group = factor(group)), 
          savePlots = TRUE, 
          suppressInds = TRUE,
          outputName = fname)

# stacked runs
temp <- runs
temp$group = "Golden"

fname = paste(config$prefix, "_stacked_runs.pdf", sep="")
fname = file.path(config$base_folder, config$outdir, fname)

plot_StackedRuns(runs = temp, savePlots = TRUE, outputName = fname)

## manhattan plots
fname = paste(config$prefix, "_manhattanRuns.pdf", sep="")
fname = file.path(config$base_folder, config$outdir, fname)

plot_manhattanRuns(
    runs = temp %>% mutate(chrom = as.character(chrom)),
    genotypeFile = pedfilePath, 
    mapFile = mapfilePath, 
    # pct_threshold = config$thr_value, 
    savePlots = FALSE, 
    outputName = fname)

## plot SNP in runs
# plot_SnpsInRuns(runs = filter(runs, chrom == 12) %>% mutate(chrom = as.character(chrom)), genotypeFile = pedfilePath, mapFile = mapfilePath)

### HRR islands
writeLines(' - calculating HRR islands')
# thr_value = 0.20
fname = paste(config$base_folder, config$repo_name, "/scripts/replacement_functions.R", sep = "")
source(fname)
fname = paste(config$base_folder, config$outdir, config$prefix, "_hrr_islands.csv", sep="")
tableRuns(runs = mutate(runs, chrom = as.character(chrom)), genotypeFile = pedfilePath, mapFile = mapfilePath, threshold = config$thr_value) %>% 
  fwrite(fname, sep=",", col.names = TRUE)

### percentage of times each SNP falls in a HRR
runs_dfx = mutate(runs, chrom = as.character(chrom))
names(runs_dfx) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

mapChrom = fread(mapfilePath)
names(mapChrom) <- c("CHR","SNP_NAME","x","POSITION")

res = NULL
for (chr in unique(mapChrom$CHR)) {
  
  print(paste("analysing chromosome", chr))
  rrnn = filter(runs_dfx, CHROMOSOME == chr)
  mmpp = filter(mapChrom, CHR == chr)
  
  temp <- snpInsideRuns(runsChrom = rrnn, mapChrom = mmpp, genotypeFile = pedfilePath)
  res = rbind.data.frame(res,temp)
}

print(head(res,20))

print("DONE!!")
