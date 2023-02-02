rm(list=ls())

suppressPackageStartupMessages({

  library(svglite)
  library(pheatmap)
  library(RColorBrewer)
  library(ggtree)
  library(ggnewscale)
  library(ggtreeExtra)
  library(ggplot2)

})

args=commandArgs(trailingOnly = TRUE)
snp_phylotree<-args[1]
summary_resistance_virulence<-args[2]

cluster_file<-read.csv(list.files(pattern = "csv"),sep='')
cluster_file$cluster<-as.character(cluster_file$cluster)
cluster_file$STRAIN_ID<-as.character(cluster_file$STRAIN_ID)

logfile<-read.delim("PDOR.log",header = F,sep='\t')
ref_name <- logfile[which(logfile$V1=="ref"),"V2"]
#tree <- read.tree("SNP_alignment.treefile")
tree <- read.tree(snp_phylotree)

add<-as.data.frame(tree$tip.label)
colnames(add)<-"STRAIN_ID"

mm<-merge(cluster_file,add,by="STRAIN_ID",all = T)

mm<-mm[which(mm$STRAIN_ID %in% tree$tip.label),]

rr<-unique(rbind(mm,cluster_file))
rownames(rr)<-rr$STRAIN_ID
rr$STRAIN_ID<-NULL

report<-read.delim(summary_resistance_virulence,sep='\t')
colnames(report)[17:18]<-c("COVERAGE","IDENTITY")

report$COVERAGE<-as.character(report$COVERAGE)
report$IDENTITY<-as.character(report$IDENTITY)

ref_name<-sub(".*\\/", "", ref_name)

report<-report[which(report$Name!="Name"),]
genes_list<-grep("Gene",report$Gene.symbol,value = T,invert = T)

report<-report[which(as.character(report$Gene.symbol) %in% genes_list),]


#new_ref_name<-paste(unique(report[which(report$X.FILE==ref_name),"SEQUENCE"]),".fna",sep="")
#new_ref_name<-unique(report[which(report$X.FILE==ref_name),"SEQUENCE"])

#report$X<-gsub(ref_name,new_ref_name,report$X.FILE)

#report[which(report$Gene.symbol=="blaKPC-3"),"Name"]

report<-report[grep("fna|fasta",report$Name),]

report$Name<-gsub(".fna","",report$Name)
report$Name<-gsub(".fasta","",report$Name)
# report$X.COVERAGE<-as.numeric(as.character(report$X.COVERAGE))
# report$X.IDENTITY<-as.numeric(as.character(report$X.IDENTITY))

report_res<-report[which(report$COVERAGE=="100.00" & report$IDENTITY=="100.00"),]
report_res<-report_res[which(report_res$Element.type=="AMR" | report_res$Element.type=="STRESS"),]
report_vir<-report[which(as.numeric(report$COVERAGE)>80 & as.numeric(report$IDENTITY)>80),]
report_vir<-report_vir[which(report_vir$Element.type=="VIRULENCE"),]

report_summary<-rbind(report_res,report_vir)

pp<-table(report_summary$Name,report_summary$Gene.symbol)


###
res_tab<-table(report_res$Name,report_res$Gene.symbol)
res_mtab<-as.matrix(res_tab)
res_mtab[res_mtab>1]<-1

res_mtab1<-matrix(as.character(res_mtab),ncol = ncol(res_mtab))
rownames(res_mtab1)<-rownames(res_mtab)
colnames(res_mtab1)<-colnames(res_mtab)

res_pres<-unique(report_res$Name)
missed_res<-unique(report[-which(report$Name %in% res_pres),"Name"])

ss<-data.frame(matrix(NA,nrow=length(missed_res),ncol=ncol(res_mtab1)))
colnames(ss)<-colnames(res_mtab1)
rownames(ss)<-missed_res
ss[is.na(ss)]<-0

res_mtab1<-rbind(res_mtab1,ss)

res_mtab1<-as.matrix(res_mtab1)

res_mtab1[res_mtab1=="1"]<-"AMR/STRESS"

###

vir_tab<-table(report_vir$Name,report_vir$Gene.symbol)
vir_mtab<-as.matrix(vir_tab)
vir_mtab[vir_mtab>1]<-1
vir_mtab1<-matrix(as.character(vir_mtab),ncol = ncol(vir_mtab))
rownames(vir_mtab1)<-rownames(vir_mtab)
colnames(vir_mtab1)<-colnames(vir_mtab)

vir_pres<-unique(report_vir$Name)
missed_vir<-unique(report[-which(report$Name %in% vir_pres),"Name"])

ss<-data.frame(matrix(NA,nrow=length(missed_vir),ncol=ncol(vir_mtab1)))
colnames(ss)<-colnames(vir_mtab1)
rownames(ss)<-missed_vir
ss[is.na(ss)]<-0

vir_mtab1<-rbind(vir_mtab1,ss)


vir_mtab1$`Gene symbol`<-NULL
vir_mtab1<-as.matrix(vir_mtab1)
vir_mtab1[vir_mtab1=="1"]<-"VIR"

summary_resvir<-cbind(vir_mtab1,res_mtab1)
summary_resvir[summary_resvir=="0"]<-NA



###

add<-as.data.frame(tree$tip.label)
colnames(add)<-"STRAIN_ID"



query_out_lab<-c()
for (i in cluster_file$cluster){
  cc<-cluster_file[which(cluster_file$cluster==i),]
  if (length(grep("DB",cc$STRAIN_ID,invert = T,value = T))>1)
    query_out_lab<-unique(c(query_out_lab,i))

}

cluster_file_query_out<-cluster_file[which(cluster_file$cluster %in% query_out_lab),]


mm<-merge(cluster_file_query_out,add,by="STRAIN_ID",all = T)
mm<-mm[which(mm$STRAIN_ID %in% tree$tip.label),]


rr<-unique(rbind(mm,cluster_file_query_out))

p<-ggtree(tree,branch.length = "none",layout = "rectangular",size=2)

rr$cluster<-ifelse(is.na(rr$cluster),"Other",rr$cluster)

rr<-rr[order(rr$cluster),]

ncols<-nrow(table(rr$cluster))

cols<-brewer.pal(n = nrow(table(rr$cluster)), name = "Dark2")

color_vec <- cols
names(color_vec) <- unique(rr$cluster)
color_vec[names(color_vec) == "Other"] <- "lightgrey"


pp <- p %<+% rr + geom_tiplab(aes(fontface=length(rr$STRAIN_ID)*0.005),alpha=1,size=length(rr$STRAIN_ID)*0.08,show.legend = FALSE)+
  theme_tree2()+theme(legend.text=element_text(size=length(rr$STRAIN_ID)*0.5),legend.title = element_text(size=length(rr$STRAIN_ID)*0.5))

nclust<-length(rr$STRAIN_ID)*0.1




res_hm <- gheatmap(pp,summary_resvir, offset = 5, width=ncol(res_mtab1)*0.3/16, font.size=nclust, colnames_position= "top",
                  colnames_angle = 90, colnames_offset_y = -0.2, hjust = 0) +
  #scale_fill_manual(values = c("0" = "white","1" = "#0F52BA"))+
  scale_y_continuous(limits=c(-1, NA))+ #guides(fill="none")+
  scale_fill_manual(breaks=c("AMR/STRESS", "VIR"),
                    values=c("#6495ED", "firebrick"), name="Category",na.translate = F)


#res_vir_hm_level <-res_hm + new_scale_fill()

rownames(rr)<-rr$STRAIN_ID
rr$STRAIN_ID<-NULL


cluster_hm_levels <-res_hm + new_scale_fill()

hm_cluster <- gheatmap(cluster_hm_levels,rr, offset = 3, width=0.02,font.size = nclust)+
  scale_fill_manual(values = color_vec,na.translate = F,name="Cluster")+
  scale_y_continuous(limits=c(-1, NA))+theme(legend.key.size = unit(3, 'cm'),legend.spacing.y = unit(5,"cm"))


ggsave(hm_cluster, filename = "annotated_tree.pdf",width = length(tree$tip.label)*1.5, height=length(tree$tip.label)*1.7 ,
       units = "cm",scale=1.2,limitsize = F)



#ggsave(hm_cluster, filename = "annotated_tree.eps",width = 70, height=70,
#      units = "cm",scale=1.5)
