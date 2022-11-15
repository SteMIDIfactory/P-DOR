rm(list=ls())

suppressPackageStartupMessages({
    
library(svglite)
library(pheatmap)
library(RColorBrewer)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(reshape2)
library(ggplot2)

})


args=commandArgs(trailingOnly = TRUE)
snp_phylotree<-args[1]

cluster_file<-read.csv(list.files(pattern = "csv"),sep='')
cluster_file$cluster<-as.character(cluster_file$cluster)
cluster_file$STRAIN_ID<-as.character(cluster_file$STRAIN_ID)

logfile<-read.delim("PDOR.log",header = F,sep='\t')
ref_name <- logfile[which(logfile$V1=="ref"),"V2"]
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

ref_name<-sub(".*\\/", "", ref_name)

new_ref_name<-paste(unique(report[which(report$X.FILE==ref_name),"SEQUENCE"]),".fna",sep="")
#new_ref_name<-unique(report[which(report$X.FILE==ref_name),"SEQUENCE"])


report$X.FILE<-gsub(ref_name,new_ref_name,report$X.FILE)



report<-report[grep("fna|fasta",report$X.FILE),]

report$X.FILE<-gsub(".fna","",report$X.FILE)
report$X.FILE<-gsub(".fasta","",report$X.FILE)

report_res<-report[which(report$X.COVERAGE=="100.00" & report$X.IDENTITY=="100.00"),]

report_res<-report_res[grep("card|resfinder|ncbi|argannot",report_res$DATABASE),]


report_vir<-report[which(as.numeric(report$X.COVERAGE)>80 & as.numeric(report$X.IDENTITY)>80),]

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

palette2 <- c("0" = "grey",
              "1" = "#960018")


add<-as.data.frame(tree$tip.label)
colnames(add)<-"STRAIN_ID"

mm<-merge(cluster_file,add,by="STRAIN_ID",all = T)
mm<-mm[which(mm$STRAIN_ID %in% tree$tip.label),]


rr<-unique(rbind(mm,cluster_file))

p<-ggtree(tree,branch.length = "none",layout = "rectangular")
rr$cluster<-ifelse(is.na(rr$cluster),"Other",rr$cluster)

rr<-rr[order(rr$cluster),]

ncols<-nrow(table(rr$cluster))
cols <- colorRampPalette(brewer.pal(8, "Dark2"))(ncols)
#cols<-brewer.pal(n = nrow(table(rr$cluster)), name = "Dark2")
color_vec <- cols
names(color_vec) <- unique(rr$cluster)
color_vec[names(color_vec) == "Other"] <- "darkgrey"



pp <- p %<+% rr + geom_tiplab(aes(color=cluster,fontface=2),alpha=1,size=length(rr$STRAIN_ID)*12/100,show.legend = FALSE)+
  geom_tippoint(aes(color=cluster),size=length(rr$STRAIN_ID)*7/100,alpha = 0)+
  theme_tree2()+theme(legend.text=element_text(size=length(rr$STRAIN_ID)*40/100),legend.title = element_text(size=length(rr$STRAIN_ID)*40/100))+
  guides(colour = guide_legend("Cluster", override.aes = list(size = length(rr$STRAIN_ID)*40/100, alpha = 1)))+
  scale_color_manual(values = color_vec)



suppressMessages(hm <- gheatmap(pp,mtab1, offset = 2, width=0.3, font.size=length(rr$STRAIN_ID)*12/100, colnames_position= "top",
               colnames_angle = 90, colnames_offset_y = -0.2, hjust = 0) +
  scale_fill_manual(values = c("0" = "white","1" = "#960018"))+
  scale_y_continuous(limits=c(-1, NA))+ guides(fill="none"))



ggsave(hm, filename = "annotated_tree.svg",width = length(tree$tip.label)*1.5, height=length(tree$tip.label)*1.5 ,
       units = "cm",scale=1.3,limitsize = F)



