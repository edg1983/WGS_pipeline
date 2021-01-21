##################################################
# Author: Edoardo Giacopuzzi
# 
# Part of the QC_bam.nf nextflow pipeline
# This R script generates plots and tables summarizing coverage and mapping data
# Expected input files in the input directory
# - flagstat (.flagstat)
# - mapping stats (.mapstats)
# - mosdepth coverage output (.global.dist.total.txt + region.dist.total.txt
##################################################

#Libraries
library(ggplot2)
library(tidyr)
library(gridExtra)
library(rmarkdown)
library(ggrepel)

#Read arguments
args = commandArgs(trailingOnly=TRUE)
Rmd_file=args[1]
assembly=args[2] #hg19 or hg38
inputdir=args[3]
ped_file=args[4]

#Standard chrs names
std_chrs <- c(paste0("chr",seq(1,22,1)), "chrX","chrY","chrM","contigs","*")

################################
###  Analyze global coverage ###
################################

#Load files and set samples names
files <- list.files(path=inputdir, pattern="global.dist.total.txt$")
cumcov <- read.table(paste0(inputdir,"/",files[1]), header = F, as.is=T)
cumcov$sample <- gsub("\\..+", "", files[1])
for (f in files[-1]) {
  mydf <- read.table(paste0(inputdir,"/",f), header = F, as.is=T)
  mydf$sample <- gsub("\\..+", "", f)
  cumcov <- rbind(cumcov, mydf)
}

# Cov tranches 
breaks <- c(1,5,10,20)

# Summary table of coverage
cov_summary <- data.frame(samples=gsub("\\..+", "", files), mean_cov=rep(0, length(files)), cov_1X = rep(0, length(files)), cov_10X=rep(0, length(files)), cov_20X=rep(0, length(files)), stringsAsFactors = F)
for (n in 1:nrow(cov_summary)) {
  s <- cov_summary$samples[n]
  mean_cov <- sum(cumcov$V3[cumcov$sample==s] * cumcov$V2[cumcov$sample==s]) / sum(cumcov$V3[cumcov$sample==s])
  newline <- c(mean_cov, cumcov$V3[cumcov$sample==s & cumcov$V2==1],cumcov$V3[cumcov$sample==s & cumcov$V2==10], cumcov$V3[cumcov$sample==s & cumcov$V2==20])
  cov_summary[n,2:5] <- newline
}
cov_summary$cov_10X <- as.numeric(cov_summary$cov_10X)
cov_summary$cov_20X <- as.numeric(cov_summary$cov_20X)
cov_summary$cov_category <- "GOOD"
cov_summary$cov_category[cov_summary$cov_10X < 0.8] <- "LOWCOV (10X < 0.8)"
cov_summary$cov_category[cov_summary$cov_20X < 0.8] <- "MIDCOV (20X < 0.8)"
cov_summary <- cov_summary[order(cov_summary$cov_10X),]

# Cumulative coverage plot
cumcov_plot <- ggplot(cumcov, aes(x=V2,y=V3,color=sample)) + geom_line() + geom_vline(xintercept = c(10,20), linetype="dashed") + geom_hline(yintercept = 0.90, linetype="dashed") + lims(x=c(0,200)) + labs(x="Coverage", y="Fraction of bases covered \u2265 X")

#Save to output
write.table(cov_summary, file = "genome.covsummary.tsv", sep="\t", row.names = F, quote=F)
ggsave(cumcov_plot, filename = "genome.cumcov.pdf", device="pdf", height=5, width=6, dpi = 150)

################################
###  Analyze region coverage ###
################################

#Load files and set samples names
files <- list.files(path=inputdir, pattern="region.dist.total.txt$")
reg_cumcov <- read.table(paste0(inputdir,"/",files[1]), header = F, as.is=T)
reg_cumcov$sample <- gsub("\\..+", "", files[1])
for (f in files[-1]) {
  mydf <- read.table(paste0(inputdir,"/",f), header = F, as.is=T)
  mydf$sample <- gsub("\\..+", "", f)
  reg_cumcov <- rbind(reg_cumcov, mydf)
}

# Cov tranches 
breaks <- c(1,5,10,20)

# Summary table of coverage
reg_cov_summary <- data.frame(samples=gsub("\\..+", "", files), mean_cov=rep(0, length(files)), cov_1X = rep(0, length(files)), cov_10X=rep(0, length(files)), cov_20X=rep(0, length(files)), stringsAsFactors = F)
for (n in 1:nrow(reg_cov_summary)) {
  s <- reg_cov_summary$samples[n]
  mean_cov <- sum(reg_cumcov$V3[reg_cumcov$sample==s] * reg_cumcov$V2[reg_cumcov$sample==s]) / sum(reg_cumcov$V3[reg_cumcov$sample==s])
  newline <- c(mean_cov, reg_cumcov$V3[reg_cumcov$sample==s & reg_cumcov$V2==1],reg_cumcov$V3[reg_cumcov$sample==s & reg_cumcov$V2==10], reg_cumcov$V3[reg_cumcov$sample==s & reg_cumcov$V2==20])
  reg_cov_summary[n,2:5] <- newline
}
reg_cov_summary$cov_10X <- as.numeric(reg_cov_summary$cov_10X)
reg_cov_summary$cov_20X <- as.numeric(reg_cov_summary$cov_20X)
reg_cov_summary$cov_category <- "GOOD"
reg_cov_summary$cov_category[reg_cov_summary$cov_10X < 0.8] <- "LOWCOV (10X < 0.8)"
reg_cov_summary$cov_category[reg_cov_summary$cov_20X < 0.8] <- "MIDCOV (20X < 0.8)"
reg_cov_summary <- reg_cov_summary[order(reg_cov_summary$cov_10X),]

# Cumulative coverage plot
reg_cumcov_plot <- ggplot(reg_cumcov, aes(x=V2,y=V3,color=sample)) + geom_line() + geom_vline(xintercept = c(10,20), linetype="dashed") + geom_hline(yintercept = 0.90, linetype="dashed") + lims(x=c(0,200)) + labs(x="Coverage", y="Fraction of bases covered \u2265 X")

#Save to output
write.table(reg_cov_summary, file = "regions.covsummary.tsv", sep="\t", row.names = F, quote=F)
ggsave(reg_cumcov_plot, filename = "regions.cumcov.pdf", device="pdf", height=5, width=6, dpi = 150)


##############################
###  Analyze mapping stats ###
##############################

#Load files and set samples names
files <- list.files(path=inputdir, pattern="mapstat$")

mapstats <- read.table(paste0(inputdir,"/",files[1]), header = F, as.is=T)
contigs <- setdiff(mapstats$V1, std_chrs)
contig_line <- c("contigs", sum(mapstats$V2[mapstats$V1 %in% contigs]), sum(mapstats$V3[mapstats$V1 %in% contigs]), sum(mapstats$V4[mapstats$V1 %in% contigs]))
mapstats <- rbind(mapstats, contig_line)
mapstats$sample <- gsub("\\..+", "", files[1])

for (f in files[-1]) {
  mydf <- read.table(paste0(inputdir,"/",f), header = F, as.is=T)
  contigs <- setdiff(mydf$V1, std_chrs)
  contig_line <- c("contigs", sum(mydf$V2[mydf$V1 %in% contigs]), sum(mydf$V3[mydf$V1 %in% contigs]), sum(mydf$V4[mydf$V1 %in% contigs]))
  mydf <- rbind(mydf, contig_line)
  mydf$sample <- gsub("\\..+", "", f)
  mapstats <- rbind(mapstats, mydf)
}
mapstats <- mapstats[mapstats$V1 %in% std_chrs,]
colnames(mapstats) <- c("CHR","LENGTH","MAPPED","UNMAPPED","SAMPLE")
mapstats$LENGTH <- as.numeric(mapstats$LENGTH)
mapstats$MAPPED <- as.numeric(mapstats$MAPPED)
mapstats$UNMAPPED <- as.numeric(mapstats$UNMAPPED)

mapstats$TOTAL <- mapstats$MAPPED + mapstats$UNMAPPED
mapstats$PCT_MAPPED <- mapstats$MAPPED / mapstats$TOTAL

#Aggregate counts per sample
mapstats_agg <- aggregate(mapstats[mapstats != "contigs",3:4],by = list(mapstats$SAMPLE[mapstats != "contigs"]), FUN=sum)
mapstats_agg$group <- "Std_Chromosomes"
mydf <- aggregate(mapstats[mapstats == "contigs",3:4],by = list(mapstats$SAMPLE[mapstats == "contigs"]), FUN=sum)
mydf$group <- "Contigs"
mapstats_agg <- rbind(mapstats_agg, mydf)
mydf <- aggregate(mapstats[,3:4],by = list(mapstats$SAMPLE), FUN=sum)
mydf$group <- "Total"
mapstats_agg <- rbind(mapstats_agg, mydf)
mapstats_agg$TOTAL <- mapstats_agg$MAPPED + mapstats_agg$UNMAPPED
mapstats_agg$PCT_MAPPED <- mapstats_agg$MAPPED / mapstats_agg$TOTAL
mapstats_agg$FLAG <- "GOOD"
mapstats_agg$FLAG[mapstats_agg$PCT_MAPPED < 0.90] <- "MIDMAP (mapped < 90%)" 
mapstats_agg$FLAG[mapstats_agg$PCT_MAPPED < 0.80] <- "LOWMAP (mapped < 80%)"

#Boxplot of mapping by chr
panel1 <- ggplot(mapstats[mapstats$CHR != "*" & mapstats$CHR != "chrM",], aes(x=factor(CHR, levels=std_chrs), y=MAPPED)) + geom_boxplot() + theme(axis.text.x = element_text(angle=45, hjust=1)) + labs(x="", y="Reads count")
panel2 <- ggplot(mapstats[mapstats$CHR != "*" & mapstats$CHR != "chrM",], aes(x=factor(CHR, levels=std_chrs), y=MAPPED/(LENGTH/1000000))) + geom_boxplot() + theme(axis.text.x = element_text(angle=45, hjust=1)) + labs(x="", y="Normalized read count (reads/Mbp)")
map_bychr <- arrangeGrob(panel1,panel2,top="Distribution of mapped reads across chromosomes")

#Save outputs
write.table(mapstats_agg, file = "mapping_summary.tsv", sep="\t", row.names = F, quote = F)
ggsave(map_bychr, filename = "mapping_bychr.pdf", height=6, width=8, dpi = 150)

###########################
###  Analyze flag stats ###
###########################

#Load files and set samples names
files <- list.files(path=inputdir, pattern="flagstat$")
flagstats <- data.frame(sample=character(), mapped=numeric(), duplicated=numeric(), paired=numeric(), stringsAsFactors = F)

for (f in files) {
  mydf <- read.table(paste0(inputdir,"/",f), header = F, as.is=T, sep="+", fill=NA)
  newline <- c(gsub("\\..+", "", f), mydf[5,1], mydf[4,1], mydf[9,1])
  flagstats[nrow(flagstats)+1,] <- newline
}
flagstats[2:4] <- sapply(flagstats[2:4],as.numeric)
flagstats$pct_dup <- flagstats$duplicated / flagstats$mapped
flagstats$pct_paired <- flagstats$paired / flagstats$mapped

#Summary table of dup and paired reads
write.table(flagstats, file = "dup_paired_pct.tsv", sep="\t", row.names = F, quote=F)

#################
### Sex Check ###
#################

#hg19 sizes from:
#/well/gel/HICF2/ref/genomes/hs37d5/hs37d5.fa.chrSizes.chrom.sizes
chrsizes_hg19 <- list()
chrsizes_hg19[["chr1"]] <- 249250621
chrsizes_hg19[["chrX"]] <- 155270560
chrsizes_hg19[["chrY"]] <- 59373566

#hg38 sizes from:
#/well/gel/HICF2/ref/genomes/hs37d5/hs37d5.fa.chrSizes.chrom.sizes
chrsizes_hg38 <- list()
chrsizes_hg38[["chr1"]] <- 248956422
chrsizes_hg38[["chrX"]] <- 156040895
chrsizes_hg38[["chrY"]] <- 57227415

if (assembly == "hg19") {
  nfactor_x <- chrsizes_hg19[["chrX"]]
  nfactor_y <- chrsizes_hg19[["chrY"]]
  nfactor_1 <- chrsizes_hg19[["chr1"]]
} else if (assembly == "hg38")  {
  nfactor_x <- chrsizes_hg38[["chrX"]]
  nfactor_y <- chrsizes_hg38[["chrY"]]
  nfactor_1 <- chrsizes_hg38[["chr1"]]
}

sexes <- read.table(ped_file, sep="\t",header=F,as.is=T)

files <- list.files(path=inputdir, pattern="mapstat$")
sample <- gsub("\\..+","",files[1],perl=T)
sex_check <- read.table(paste0(inputdir,"/",files[1]), sep="\t",header=F,as.is=T)
sex_check$sample <- sample

for (f in files[-1]) {
  #message(paste0("Processing ", f))
  sample <- gsub("\\..+","",f,perl=T)
  mydf <- read.table(paste0(inputdir,"/",f), sep="\t",header=F,as.is=T)
  mydf$sample <- sample
  sex_check <- rbind(sex_check,mydf)
}

sex_check$V1 <- gsub("chr","",sex_check$V1)
used_chrs <- c("1","X","Y")
sex_check <- sex_check[sex_check$V1 %in% used_chrs, ]
sex_check <- spread(sex_check[,c(1,3,5)], V1,V3)

sex_check$norm_1 <- sex_check$`1` / nfactor_1
sex_check$norm_X <- sex_check$X / nfactor_x
sex_check$norm_Y <- sex_check$Y / nfactor_y
sex_check$Xto1_ratio <- sex_check$norm_X/sex_check$norm_1
sex_check$Yto1_ratio <- sex_check$norm_Y/sex_check$norm_1

sex_check <- merge(sex_check,sexes[,c(2,5)], by.x="sample", by.y="V2")
colnames(sex_check)[10] <- "Sex"

sex_check$imputed_Sex <- "0"
sex_check$imputed_Sex[sex_check$Yto1_ratio > 0.3 & sex_check$Xto1_ratio < 0.7] <- "1"
sex_check$imputed_Sex[sex_check$Yto1_ratio < 0.3 & sex_check$Xto1_ratio > 0.7] <- "2"

sex_check$checkSex_result <- "CONCORDANT"
sex_check$checkSex_result[sex_check$Sex != sex_check$imputed_Sex & sex_check$Sex !=0 & sex_check$imputed_Sex != 0] <- "DISCORDANT"
sex_check$checkSex_result[sex_check$imputed_Sex != 0 & sex_check$Sex == 0 ] <- "IMPUTED"
sex_check$checkSex_result[sex_check$imputed_Sex == 0] <- "UNCERTAIN"

write.table(sex_check[,c(1,8:12)], file = "Sex_check_results.tsv", sep="\t", row.names=F, quote=F)

#####################################
### Save R data and generate html ###
#####################################

save.image("QC_bam.Rdata")
#rmarkdown::render("QC_bam.Rmd", "html_document", "QC_bam.html")
