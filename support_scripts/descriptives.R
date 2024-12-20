library("ggplot2")
library("tidyverse")
library("data.table")

basefolder = "/home/filippo/Documents/ciampolini/golden_retriever"
infile_miss_before = "results/golden.imiss"
infile_freq_after = "results/golden_filtered.frq"
infile_ld = "results/ld/golden_filtered.ld"
infile_ld_chr = "results/ld/golden_filtered_38.ld"

## before filtering
fname = file.path(basefolder, infile_miss_before)
imiss = fread(fname)

imiss |>
  summarise(avg = mean(F_MISS), std = sd(F_MISS))

## after filtering
## missing rate is 0% given the current sample size (N = 18) and the chosen per-SNP threshold
## of max 5% missing rate allowed

fname = file.path(basefolder, infile_freq_after)
frq = fread(fname)

###################################################
## min n. of SNP in ROH (from Mastrangelo et al. 2018)
q = frq$MAF
p = 1-q
1-mean(2*p*q)
mean(p**2)
mean(q**2)

alpha = 0.01
n = nrow(imiss)
lgen = nrow(frq)

s = log(alpha/(n*lgen)) / (log(1-(1-mean(2*p*q))))
###################################################

hist(frq$MAF, main = "MAF", col = "darkgreen")

## LD

# library("devtools")
# devtools::install_github('Rong-Zh/GWLD/GWLD-R')
library("GWLD")

# data(duck)
# data <- duck$SNP
# SNP <- data$genotype
# Info <- data$info
# 
# SNP <- codegeno(SNP, sep="/")
# 
# result <- GWLD(SNP, method = "r^2", cores = 1)
# 
# p <- HeatMap(data = SNP, method = "RMI", SnpPosition = Info$POS, SnpName = Info$ID, cores = 1, color = "YellowTored", showLDvalues = FALSE)

library("gaston")
# LD.plot(LD = result)

fname = file.path(basefolder, infile_ld_chr)
ldres = fread(fname)

x <- ldres[6000:7000,] 

temp <- x |>
  select(SNP_A, SNP_B, R2) |>
  tidyr::spread(key = "SNP_B", value = "R2", fill = 0) |>
  select(-SNP_A) |>
  as.matrix()

# LD.plot(result)
LD.plot(temp)

###########
## LD decay
###########
## Golden
D <- data.frame("dbin"=NULL,"media"=NULL,"dist"=NULL,"n"=NULL,"breed"=NULL)

fname = file.path(basefolder, infile_ld)
linkd <- fread(fname, header = TRUE)
linkd$distance <- abs(linkd$BP_A-linkd$BP_B)
  
linkd$dbin = findInterval(linkd$distance,seq(0,max(linkd$distance),50*10^3))

linkd |>
  select(R2, distance) |>
  summarise(across(where(is.numeric), list(avg = mean, std = sd, med = median))) |>
  fwrite("results/ld/ld.csv", sep = ",", col.names = TRUE)

medias <- linkd %>%
  group_by(dbin) %>%
  summarize(media=median(R2),dist=mean(distance),n=n()) %>%
  filter(n>10)

breed = "golden"
medias$breed <- rep(breed,nrow(medias))
print(medias)
D <- rbind.data.frame(D,medias)

p <- ggplot()
p <- p + geom_line(data = D, aes(dist/1000, media, colour=breed), linewidth = 1.5)
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom", legend.text=element_text(size=5))
p <- p + xlab("Distance in kpbs") + ylab("LD as r2")
p <- p + ggtitle("Median LD as a function of distance")
p <- p + scale_x_continuous(breaks=as.vector(round(quantile(D$dist/1000,c(0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.999)))))
# p <- p + scale_colour_manual(values = heat.colors(33))
# p <- p + scale_colour_manual(values = as.character(farben$V1))
print(p)

ggsave(filename = "results/ld/LD_decay_golden.png", plot = p, device = "png")

### relative to other breeds

avgld = fread("/home/filippo/Documents/ciampolini/razzeCanine/ld/avgLD.txt")
avgld = filter(avgld, breed != "GRe")
avgld <- bind_rows(D,avgld)

D_GO <- avgld %>%
  filter(breed=="golden")

D_other <- avgld %>%
  filter(breed != "golden")

D_other$breed <- toupper(D_other$breed)

farben <- read.table("~/Documents/ciampolini/razzeCanine/Ne/colors.txt", header = FALSE)

p <- ggplot()
p <- p + geom_line(data = D_other, aes(dist/1000, media, colour=breed))
p <- p + geom_line(data = D_GO, aes(dist/1000, media, colour=breed), size=1.5)
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom", legend.text=element_text(size=5))
p <- p + xlab("Distance in kpbs") + ylab("LD as r2")
p <- p + ggtitle("Median LD as a function of distance")
p <- p + scale_x_continuous(breaks=as.vector(round(quantile(D$dist/1000,c(0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.999)))))
# p <- p + scale_colour_manual(values = heat.colors(33))
# p <- p + scale_colour_manual(values = as.character(farben$V1))
print(p)

ggsave(filename = "results/ld/LD_decay.png", plot = p, device = "png", width = 8, height = 10)


