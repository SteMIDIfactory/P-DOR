rm(list=ls())
library(ape)
library(openxlsx)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)


args=commandArgs(trailingOnly = TRUE)
snp_alignment <- args[1]
snp_threshold <- args[2]

#file_folder <- list.files()
#file_folder<-file_folder[grep("Results",file_folder)]
###
seq<-read.dna(snp_alignment, format="fasta")
dist_matrix<-as.matrix(dist.dna(seq, "N"))
###
query_genomes<-as.vector(grep("DB_",rownames(dist_matrix),value = TRUE, invert = TRUE))
# background_genomes<-as.vector(grep("DB_",rownames(dist_matrix),value = TRUE))
ref_genome<-as.vector(grep("REF_",rownames(dist_matrix),value = TRUE))
query_genomes<-query_genomes[-which(query_genomes %in% ref_genome)]
query_genomes<-dist_matrix[query_genomes,query_genomes] #relativa ai genomi di outbreak di studio
# background<-dist_matrix[background_genomes,background_genomes] #relativa ai genomi di outbreak di studio


svg(filename = "SNPs_frequency.svg",
    width = 10, height = 10, pointsize = 8)
snp_plot<-hist(dist_matrix,breaks=200)
snp_plot<-hist(dist_matrix,breaks=200, xlim = range(0, max(dist_matrix)+max(dist_matrix)*0.2),
               ylim=range(0, max(snp_plot$counts)+max(snp_plot$breaks)*0.05),col = c("darkgreen"),
               xlab = "pairwise SNPs distance")
abline(v = args[2], col="red", lwd=3, lty=2)
dev.off()

### INFLECTION POINT
d1 <- density(dist_matrix)
DeltaY = diff(d1$y)
Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1

svg(filename = "inflection_points.svg",
    width = 10, height = 10, pointsize = 8)
plot(d1, xlab="", ylab="", main="")
points(d1$x[Turns], d1$y[Turns], pch=16, col="red")
dev.off()


dist_matrix_query_vs_all<-dist_matrix[rownames(query_genomes),]

###
svg(filename = "SNP_heatmap_query_vs_all.svg",
   width = 10, height = 10, pointsize = 8)
query_vs_background<-pheatmap(dist_matrix_query_vs_all,color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100))
dev.off()
