rm(list=ls())

library(svglite)
library(pheatmap)
library(RColorBrewer)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(reshape2)
library(ggplot2)
args=commandArgs(trailingOnly = TRUE)
snp_phylotree<-args[1]

cluster_file<-read.csv(list.files(pattern = "csv"),sep='')
cluster_file$cluster<-as.character(cluster_file$cluster)
cluster_file$STRAIN_ID<-as.character(cluster_file$STRAIN_ID)

logfile<-read.delim("PDOR.log",header = F,sep='\t')
ref_name <- logfile[which(logfile$V1=="ref"),"V2"]
ref_name
#tree <- read.tree("RAxML_bipartitionsBranchLabels.coreSNPs_Phylo.nwk")
tree <- read.tree(snp_phylotree)

add<-as.data.frame(tree$tip.label)
colnames(add)<-"STRAIN_ID"

mm<-merge(cluster_file,add,by="STRAIN_ID",all = T)

mm<-mm[which(mm$STRAIN_ID %in% tree$tip.label),]

rr<-unique(rbind(mm,cluster_file))
rownames(rr)<-rr$STRAIN_ID
rr$STRAIN_ID<-NULL


report<-read.delim("summary_resistance_virulence",sep='\t')

report$X.COVERAGE<-as.character(report$X.COVERAGE)
report$X.IDENTITY<-as.character(report$X.IDENTITY)
new_ref_name<-paste(unique(report[which(report$X.FILE==ref_name),"SEQUENCE"]),".fna",sep="")

report$X.FILE<-gsub(ref_name,new_ref_name,report$X.FILE)

report<-report[grep("fna|fasta",report$X.FILE),]

report$X.FILE<-gsub(".fna","",report$X.FILE)
report$X.FILE<-gsub(".fasta","",report$X.FILE)

report_res<-report[which(report$X.COVERAGE=="100.00" & report$X.IDENTITY=="100.00"),]
report_res<-report_res[grep("card|resfinder|ncbi|argannot",report_res$DATABASE),]

report_vir<-report[which(as.numeric(report$X.COVERAGE)>90 & as.numeric(report$X.IDENTITY)>90),]
report_vir<-report_vir[grep("ecoli_vf|vfdb",report_vir$DATABASE),]

report_summary<-rbind(report_res,report_vir)
tab<-table(report_summary$X.FILE,report_summary$GENE)
# tab[tab==0]<-"absence"
# tab[tab>=1]<-"presence"
mtab<-as.matrix(tab)
mtab[mtab>1]<-1
colnames(mtab)<-gsub("Bla","",colnames(mtab))
colnames(mtab)<-gsub("[(]","",colnames(mtab))
colnames(mtab)<-gsub("[)]","",colnames(mtab))

mtab1<-matrix(as.character(mtab),ncol = ncol(mtab))
rownames(mtab1)<-rownames(mtab)
colnames(mtab1)<-colnames(mtab)

palette2 <- c("0" = "grey95",
              "1" = "#960018")


add<-as.data.frame(tree$tip.label)
colnames(add)<-"STRAIN_ID"

mm<-merge(cluster_file,add,by="STRAIN_ID",all = T)
mm<-mm[which(mm$STRAIN_ID %in% tree$tip.label),]

rr<-unique(rbind(mm,cluster_file))

p<-ggtree(tree,branch.length = "none",layout = "rectangular")
rr$cluster<-ifelse(is.na(rr$cluster),"Other",rr$cluster)

rr<-rr[order(rr$cluster),]

cols<-brewer.pal(n = nrow(table(rr$cluster)), name = "Paired")
color_vec <- cols
names(color_vec) <- unique(rr$cluster)
color_vec[names(color_vec) == "Other"] <- "darkgrey"



pp <- p %<+% rr + geom_tiplab(aes(color=cluster,fontface=2),size=5,show.legend = FALSE)+
  geom_tippoint(aes(color=cluster),size=5,alpha = 0)+#scale_color_discrete(breaks = levels(factor(rr$cluster)))+
  theme_tree2()+theme(legend.text=element_text(size=20),legend.title = element_text(size=20))+
  guides(colour = guide_legend("Cluster", override.aes = list(size = 20, alpha = 1)))+
  scale_color_manual(values = color_vec)#+geom_fruit(geom=geom_point,mapping=aes(color=cluster),size=3)


hm <- gheatmap(pp,mtab1, offset = 2, width=0.2, font.size=5, colnames_position= "top",
         colnames_angle = 90, colnames_offset_y = -0.2, hjust = 0) +
         scale_fill_manual(values = c("0" = "grey95","1" = "#960018"))+
         scale_y_continuous(limits=c(-1, NA))+ guides(fill="none")


                                                                   
ggsave(hm, filename = "annotated_tree.svg",width = length(tree$tip.label)*1.5, height=length(tree$tip.label)*2.5 ,
        units = "cm",scale=2,limitsize = F)

