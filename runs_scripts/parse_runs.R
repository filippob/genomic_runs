# devtools::install_github("bioinformatics-ptp/detectRUNS/detectRUNS")

library("ggplot2")
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
    thr_value = 0.5, ##threshold for ROH/HRR islands
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
print(paste("HRR/ROH threshold:", config$thr_value))

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
  summarise(N=n(), n_ind = length(unique(id)), avg_n = N/n_ind, avg_length = mean(lengthBps)) %>%
  fwrite(fname, sep = ",", col.names = TRUE)

## manhattan plots
fname = paste(config$prefix, "_manhattanRuns", sep="")
fname = file.path(config$base_folder, config$outdir, fname)

plot_manhattanRuns(
    runs = runs %>% mutate(chrom = as.character(chrom)),
    genotypeFile = pedfilePath, 
    mapFile = mapfilePath, 
    pct_threshold = config$thr_value,
    savePlots = TRUE, 
    outputName = fname)


## plot runs per samples
writeLines(' - plotting individual runs')

temp <- filter(runs, chrom == 20)

fname = paste(config$prefix, "_plot_runs_chr20", sep="")
fname = file.path(config$base_folder, config$outdir, fname)
plot_Runs(runs = mutate(temp, group = factor(group)), 
          savePlots = TRUE, 
          suppressInds = TRUE,
          outputName = fname)

# stacked runs

fname = paste(config$prefix, "_stacked_runs_chr20", sep="")
fname = file.path(config$base_folder, config$outdir, fname)

plot_StackedRuns(runs = temp, savePlots = TRUE, outputName = fname)

## plot SNP in runs
# plot_SnpsInRuns(runs = filter(runs, chrom == 12) %>% mutate(chrom = as.character(chrom)), genotypeFile = pedfilePath, mapFile = mapfilePath)

### HRR islands
writeLines(' - calculating HRR islands')
# thr_value = 0.20
fname = file.path(config$base_folder, config$repo_name, "support_scripts/replacement_functions.R")
source(fname)

fname = paste(config$prefix, "_hrr_islands.csv", sep="")
fname = file.path(config$base_folder, config$outdir, fname)
detectRUNS::tableRuns(runs = mutate(runs, chrom = as.character(chrom)), 
                      mapFile = mapfilePath, 
                      genotypeFile = pedfilePath, 
                      threshold = config$thr_value) |>
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
  
  temp <- detectRUNS:::snpInsideRuns(runsChrom = rrnn, mapChrom = mmpp, genotypeFile = pedfilePath)
  res = rbind.data.frame(res,temp)
}

print(head(res,20))

### Genomic inbreeding

snps <- fread(mapfilePath)
names(snps) <- c("chrom", "snp", "cM", "bps")

genome_size = snps |>
  group_by(chrom) |>
  summarise(length = max(bps)) |>
  summarise(tot_length = sum(length)) |>
  pull(tot_length)

froh <- runs |>
  group_by(id) |>
  summarise(tot_roh = sum(lengthBps), f_roh = tot_roh/genome_size)

summary(froh$f_roh)
sd(froh$f_roh)

library("ggalt")
 p<- ggplot(froh, aes(x = reorder(id,f_roh), y = f_roh)) + 
  geom_lollipop(aes(color = id), size = 1.5) + 
  guides(color="none") + coord_flip() +
  theme(axis.text.y = element_blank())

print(p)

fname = paste(config$prefix, "_FROH.png", sep="")
fname = file.path(config$base_folder, config$outdir, fname)

ggsave(filename = fname, plot = p, device = "png")

print("DONE!!")
