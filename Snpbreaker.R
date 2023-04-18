
rm(list=ls())

suppressPackageStartupMessages({

library(ape)
library(reshape2)
library(ggplot2)
library(igraph)
library(pheatmap)
library(svglite)
library(RColorBrewer)

})

options(warn=-1)

args=commandArgs(trailingOnly = TRUE)

logfile<-read.delim("PDOR.log",header = F,sep='\t')



snp_alignment <- args[1]
snp_threshold <- logfile[which(logfile$V1=="snp_threshold"),"V2"]

###

seq<-read.dna(snp_alignment, format="fasta")
dist_matrix<-as.matrix(dist.dna(seq, "N"))

svg(filename = "SNP_heatmap.svg",
    width = ncol(dist_matrix)*0.7, height = ncol(dist_matrix)*0.7, pointsize = 8)
query_vs_background<-pheatmap(dist_matrix,color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100),fontsize=ncol(dist_matrix)*0.3,cellwidth=ncol(dist_matrix)*0.5,
	cellheight=ncol(dist_matrix)*0.5)
invisible (dev.off())


write.csv(dist_matrix,"snp_distance_matrix.tsv",sep='\t')


if (snp_threshold !="infl"){
  snp_threshold<-as.integer(as.character(snp_threshold))
  
  svg(filename = sprintf("SNP_frequency_manual_threshold_%s_snps.svg",snp_threshold),
      width = 10, height = 10, pointsize = 8)
  snp_plot<-hist(dist_matrix,breaks=2000)
  snp_plot<-hist(dist_matrix,breaks=2000, xlim = range(0, max(dist_matrix)+max(dist_matrix)*0.2),
                 ylim=range(0, max(snp_plot$counts)+max(snp_plot$breaks)*0.05),col = c("darkgreen"),
                 xlab = "pairwise SNPs distance")
  abline(v = snp_threshold, col="red", lwd=3, lty=2)
  invisible (dev.off()) 
  
  
  
  mm<-melt(dist_matrix)
  mm<-mm[which(mm$value<snp_threshold & mm$Var1!=mm$Var2),]
  links<-mm[,c(1,2)]
  
  svg(filename = sprintf("SNP_clusters_manual_threshold_%s_snps.svg",snp_threshold),
      width = 10, height = 10, pointsize = 8)
  gr<-graph_from_data_frame(links, directed = FALSE)
  plot(gr)
  invisible (dev.off())
  out_cluster<-as.data.frame(components(gr)$membership)
  colnames(out_cluster)<-"cluster"
  out_cluster$cluster<-paste("C",out_cluster$cluster,sep="")
  out_cluster$ID<-rownames(out_cluster)
  rownames(out_cluster)<-1:nrow(out_cluster)
  colnames(out_cluster)[2]<-"STRAIN_ID"
  write.table(out_cluster,sprintf("clusters_manual_threshold_%s_snps.csv",snp_threshold),row.names = F,quote = F,sep='\t')
  
} else {
  
  
  
  
  ### calculate inflection point
  d1 <- density(dist_matrix)
  DeltaY = diff(d1$y)
  Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0)-1
  min_flex<-Turns[1]
  
  
  svg(filename = sprintf("SNP_frequency_calculated_threshold_%s_snps.svg",min_flex),
      width = 10, height = 10, pointsize = 8)
  snp_plot<-hist(dist_matrix,breaks=200)
  snp_plot<-hist(dist_matrix,breaks=200, xlim = range(0, max(dist_matrix)+max(dist_matrix)*0.2),
                 ylim=range(0, max(snp_plot$counts)+max(snp_plot$breaks)*0.05),col = c("darkgreen"),
                 xlab = "pairwise SNPs distance")
  abline(v = min_flex, col="blue", lwd=3, lty=2)
  invisible (dev.off())
  
  
  
  
  mm_man_thresh<-melt(dist_matrix)
  mm_man_thresh<-mm_man_thresh[which(mm_man_thresh$value<min_flex & mm_man_thresh$Var1!=mm_man_thresh$Var2),]
  links<-mm_man_thresh[,c(1,2)]
  
  
  svg(filename = sprintf("SNP_clusters_calculated_threshold_%s_snps.svg",min_flex),
      width = 10, height = 10, pointsize = 8)
  gr<-graph_from_data_frame(links, directed = FALSE)
  out_cluster<-as.data.frame(components(gr)$membership)
  colnames(out_cluster)<-"cluster"
  out_cluster$cluster<-paste("C",out_cluster$cluster,sep="")
  out_cluster$ID<-rownames(out_cluster)
  rownames(out_cluster)<-1:nrow(out_cluster)
  colnames(out_cluster)[2]<-"STRAIN_ID"
  write.table(out_cluster,sprintf("clusters_calculated_threshold_%s_snps.csv",min_flex),row.names = F,quote = F,sep='\t')
  plot(gr)
  invisible(dev.off())
  
  
} 
