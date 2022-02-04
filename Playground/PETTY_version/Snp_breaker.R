rm(list=ls())
library(ape)
library(openxlsx)
library(ggplot2)

#setwd("/home/ghepard/Miniskynet/")
###
#args=commandArgs(trailingOnly = TRUE)
#setwd(getwd())
#snp_alignment <- args[1]
###
#system(sprintf("mkdir -p results_%s",format(Sys.time(), "%a_%b_%d_%X_%Y")))
file_folder <- list.files()
file_folder<-file_folder[grep("Results",file_folder)]
###
seq<-read.dna("SNP_alignment.core.fasta", format="fasta")

dist_matrix<-as.matrix(dist.dna(seq, "N"))
###
outbreak_genomes<-as.vector(grep("DB_",rownames(dist_matrix),value = TRUE, invert = TRUE))
external_genomes<-as.vector(grep("DB_",rownames(dist_matrix),value = TRUE))
ref_genomes<-as.vector(grep("REF_",rownames(dist_matrix),value = TRUE))
outbreak_genomes<-outbreak_genomes[- which(outbreak_genomes %in% ref_genomes)]
out_snp_tab<-dist_matrix[outbreak_genomes,outbreak_genomes] #relativa ai genomi di outbreak di studio
ext_snp_tab<-dist_matrix[external_genomes,external_genomes] #relativa ai genomi di outbreak di studio
###
pdf("all.pdf", height = 12, width = 15, pointsize = 12)
all<-hist(dist_matrix,breaks = 200)
all<-hist(dist_matrix,breaks=200, xlim = range(0, max(dist_matrix)+max(dist_matrix)*0.2),
         ylim=range(0, max(all$counts)+max(all$breaks)*0.25))
print(all)
dev.off()
###
cairo_ps("patric.eps", height = 12, width = 15, pointsize = 12, fallback_resolution = 300)
patric<-hist(ext_snp_tab,breaks = 200)
patric<-hist(ext_snp_tab,breaks=200, xlim = range(0, max(ext_snp_tab)+max(ext_snp_tab)*0.2),
          ylim=range(0, max(patric$counts)+max(patric$breaks)*0.25))
print(patric)
dev.off()
###
cairo_ps("outbreak.eps", height = 12, width = 15, pointsize = 12, fallback_resolution = 300)
outbreak<-hist(out_snp_tab,breaks=200)
outbreak<-hist(out_snp_tab,breaks=200, xlim = range(0, max(out_snp_tab)+max(out_snp_tab)*0.2),
             ylim=range(0, max(outbreak$counts)+max(outbreak$breaks)*0.25))
print(outbreak)
dev.off()
###

#---------------------------------------------------------------------------------------------# OVERLAP

library(reshape2)
out<-melt(out_snp_tab)
out$type<-"outbreak"
ext<-melt(ext_snp_tab)
ext$type<-"patric"
ss<-rbind(out,ext)

overlap<-ggplot(ss,aes(value))+geom_histogram(data=subset(ss,type == 'outbreak'),binwidth = 5,fill = "red",alpha=0.5)+
  geom_histogram(data=subset(ss,type == 'patric'),fill = "royalblue",alpha=0.5)+theme_light()+
  geom_vline(xintercept=20, linetype="dashed",color = "red", size=0.5)

ggsave(overlap, filename = "overlap_plot.png",device = "png",units = "cm",width = 30,height = 20)


#----------------------------------------------------------------------------------------------# HEATMAP
library(pheatmap)
library(RColorBrewer)
library(viridis)

snp_dist_out_patr<-dist_matrix[outbreak_genomes,]
###
heat_out<-pheatmap(out_snp_tab,color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100))
cairo_ps(filename = "SNP_heatmap_outbreak.eps",
           width = 10, height = 10, pointsize = 8,
         fallback_resolution = 300)
print(heat_out)
dev.off()
###
heat_out_vs_patr<-pheatmap(t(snp_dist_out_patr),color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100))
cairo_ps(filename = "SNP_heatmap_outbreak_vs_targetDB.eps",
         width = 10, height = 10, pointsize = 8,
         fallback_resolution = 300)
print(heat_out_vs_patr)
dev.off()
#----------------------------------------------------------------------------------------------# TABS
write.xlsx(as.data.frame(dist_matrix),"all_snp_distance_tab.xlsx")
write.xlsx(as.data.frame(ext_snp_tab),"patric_snp_distance_tab.xlsx")
write.xlsx(as.data.frame(out_snp_tab),"outbreak_snp_distance_tab.xlsx")
###

#system(sprintf("mkdir -p %s/RAxML",file_folder))
#system(sprintf("mkdir -p %s/SNPs_distance",file_folder))
#system(sprintf("mkdir -p %s/SNPs_distance/tables",file_folder))
#system(sprintf("mkdir -p %s/tmp",file_folder))
#system(sprintf("mkdir -p %s/Fasta_alignments",file_folder))

#system(sprintf("mv *list* %s/tmp",file_folder))
#system(sprintf("mv query_list %s/tmp",file_folder))
#system(sprintf("mv *fasta* %s/Fasta_alignments",file_folder))
#system(sprintf("mv entry_mash_selection %s/tmp",file_folder))
#system(sprintf("mv *nwk %s/RAxML",file_folder))
#system(sprintf("mv *xlsx %s/SNPs_distance/tables",file_folder))
#system(sprintf("mv *eps %s/SNPs_distance",file_folder))
#system(sprintf("mv *png %s/SNPs_distance",file_folder))

