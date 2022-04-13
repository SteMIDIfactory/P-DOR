rm(list=ls())
library(ape)
library(openxlsx)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)

#setwd("/home/ghepard/Miniskynet/")

args=commandArgs(trailingOnly = TRUE)
snp_alignment <- args[1]
snp_threshold <- args[2]
### se soglia specificata Rscript snpbreaker snpalign threshold
### altrimenti calcolare flesso

#system(sprintf("mkdir -p results_%s",format(Sys.time(), "%a_%b_%d_%X_%Y")))
setwd("/home/ghepard/PDOR/snp_distance_matrix/")
file_folder <- list.files()
file_folder<-file_folder[grep("Results",file_folder)]
###
seq<-read.dna("SNP_alignment.core.fasta", format="fasta")
dist_matrix<-as.matrix(dist.dna(seq, "N"))
###
query_genomes<-as.vector(grep("DB_",rownames(dist_matrix),value = TRUE, invert = TRUE))
background_genomes<-as.vector(grep("DB_",rownames(dist_matrix),value = TRUE))
ref_genomes<-as.vector(grep("REF_",rownames(dist_matrix),value = TRUE))
query_genomes<-query_genomes[-which(query_genomes %in% ref_genomes)]
query_genomes<-dist_matrix[query_genomes,query_genomes] #relativa ai genomi di outbreak di studio
background<-dist_matrix[background_genomes,background_genomes] #relativa ai genomi di outbreak di studio


snp_plot<-hist(dist_matrix,breaks=200)
snp_plot<-hist(dist_matrix,breaks=200, xlim = range(0, max(dist_matrix)+max(dist_matrix)*0.2),
               ylim=range(0, max(outbreak$counts)+max(outbreak$breaks)*0.05),col = c("darkgreen"),
               xlab = "pairwise SNPs distance")
abline(v = args[2], col="red", lwd=3, lty=2)

svg(filename = "SNPs_frequency.svg",
    width = 10, height = 10, pointsize = 8)
dev.off()

### INFLECTION POINT
d1 <- density(dist_matrix)
plot(d1)
DeltaY = diff(d1$y)
Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1

svg(filename = "inflection_points.svg",
    width = 10, height = 10, pointsize = 8)
plot(d1, xlab="", ylab="", main="")
points(d1$x[Turns], d1$y[Turns], pch=16, col="red",size=10)
dev.off()

####
# query<-melt(query_genomes)
# query$type<-"outbreak"
# back<-melt(background)
# back$type<-"patric"
# snp_distance_matrix<-rbind(query,back)

dist_matrix<-dist_matrix[rownames(query_genomes),]
query_vs_background<-pheatmap(dist_matrix,color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100))
svg(filename = "SNP_heatmap_all_vs_all.svg",
   width = 10, height = 10, pointsize = 8)
print(allstrain)
dev.off()



