rm(list=ls())

library(svglite)
library(ggtree)
library(reshape2)
library(ggplot2)
args=commandArgs(trailingOnly = TRUE)
snp_phylotree<-args[1]

cluster_file<-read.csv(list.files(pattern = "csv"),sep='')
cluster_file$cluster<-as.character(cluster_file$cluster)
cluster_file$STRAIN_ID<-as.character(cluster_file$STRAIN_ID)

tree <- read.tree(snp_phylotree)

add<-as.data.frame(tree$tip.label)
colnames(add)<-"STRAIN_ID"

mm<-merge(cluster_file,add,by="STRAIN_ID",all = T)

mm<-mm[which(mm$STRAIN_ID %in% tree$tip.label),]

rr<-unique(rbind(mm,cluster_file))

p<-ggtree(tree,branch.length = "none",layout = "circular")

annotated_tree<-p %<+% rr+geom_tiplab(aes(fill=cluster,color=cluster),show.legend = FALSE,alpha=1,size=7,linesize=.5)+
 theme_tree2()+theme(legend.text=element_text(size=20),legend.title =element_text(size=25))+ geom_polygon(aes(fill = cluster, x = 0, y = 0)) +
        scale_fill_discrete(na.translate = F)



ggsave(annotated_tree, filename = "annotated_tree.svg",width = length(tree$tip.label)*1.5, height=length(tree$tip.label) , units = "cm")




