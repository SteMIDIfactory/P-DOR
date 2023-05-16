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
# summary_resistance_virulence<-args[2]


tree <- read.tree(snp_phylotree)
#tree <- read.tree("SNP_alignment.treefile")


add<-as.data.frame(tree$tip.label)
colnames(add)<-"STRAIN_ID"

####################################### CHECK CLUSTER FILE

check_cluster_file<-list.files(pattern = "csv")



if (length(check_cluster_file)>0) {
  
  cluster_file<-read.csv(list.files(pattern = "csv"),sep='')
  cluster_file$cluster<-as.character(cluster_file$cluster)
  cluster_file$STRAIN_ID<-as.character(cluster_file$STRAIN_ID)
  
  cat("\n\nPossible outbreak detected!\n\nGenerating annotated phylogenetic tree...\n\n")
  
} else {
  
  cluster_file<-data.frame(matrix(ncol=2))
  colnames(cluster_file)<-c("STRAIN_ID","cluster")
  
}

#############################################  

mm<-merge(cluster_file,add,by="STRAIN_ID",all = T)

mm<-mm[which(mm$STRAIN_ID %in% tree$tip.label),]

rr<-unique(rbind(mm,cluster_file))
rr<-rr[-which(is.na(rr$STRAIN_ID)),]

rownames(rr)<-rr$STRAIN_ID
rr$STRAIN_ID<-NULL

add<-as.data.frame(tree$tip.label)
colnames(add)<-"STRAIN_ID"

query_out_lab<-c()
for (i in cluster_file$cluster){
  cc<-cluster_file[which(cluster_file$cluster==i),]
  if (length(grep("BD_",cc$STRAIN_ID,invert = T,value = T))>=1)
    query_out_lab<-unique(c(query_out_lab,i))
  
}

cluster_file_query_out<-cluster_file[which(cluster_file$cluster %in% query_out_lab),]

mm<-merge(cluster_file_query_out,add,by="STRAIN_ID",all = T)
mm<-mm[which(mm$STRAIN_ID %in% tree$tip.label),]

rr<-unique(rbind(mm,cluster_file_query_out))

p<-ggtree(tree,branch.length = "none",layout = "rectangular",size=length(rr$STRAIN_ID)*0.025)

rr<-rr[order(rr$cluster),]

ncols<-nrow(table(rr$cluster))

rr$cluster<-ifelse(is.na(rr$cluster),"Sporadic",rr$cluster)

cols<-suppressWarnings(brewer.pal(n = nrow(table(rr$cluster)), name = "Dark2"))

if (nrow(table(rr$cluster))==1){
  cols<-cols[1]
}

color_vec <- cols
names(color_vec) <- unique(rr$cluster)
color_vec[names(color_vec) == "Sporadic"] <- "lightgrey"

labsize<-length(rr$STRAIN_ID)*0.15



########### NO CLUSTER NO RESISTANCE
if (length(unique(rr$cluster))==1 && is.na(unique(rr$cluster))){
  
  pp <- p %<+% rr + geom_tiplab(aes(fontface=length(rr$STRAIN_ID)*0.005),alpha=1,size=labsize*0.75,offset = 0.1,show.legend = FALSE)+
    theme_tree2()+theme(legend.text=element_text(size=length(rr$STRAIN_ID)*0.5),legend.title = element_text(size=length(rr$STRAIN_ID)*0.5))
  
  nclust<-length(rr$STRAIN_ID)*0.1
  
  rownames(rr)<-rr$STRAIN_ID
  rr$STRAIN_ID<-NULL
  
  
  colnames(rr)<-" "
  
  
  suppressMessages(
    
    hm_cluster <- gheatmap(pp,rr, offset = 2, width=0.02,font.size = labsize)+
      scale_fill_manual(values = color_vec,na.translate = F,name="Cluster")+
      scale_y_continuous(limits=c(-1, NA))+theme(legend.key.size = unit(nclust*0.5, 'cm'),
                                                 legend.spacing.y = unit(nclust*0.5,"cm"))
    
  )
  
  ggsave(hm_cluster, filename = "annotated_tree.svg",width = length(tree$tip.label)*1.5, height=length(tree$tip.label)*1.7 ,
         units = "cm",scale=1.2,limitsize = F)
  
  
  ########### CLUSTER YES  RESISTANCE YES
  
} else if (file.exists("summary_resistance_virulence.txt") && length(unique(rr$cluster)>=1)) {
  
  report<-read.delim("summary_resistance_virulence.txt",sep='\t')
  #kpc<-grep("bla",report$Gene.symbol,value = T)
  #report<-report[which(report$Gene.symbol %in% kpc),]
  #report$Name<-gsub("BD_","",report$Name)
  
  refname<-unique(report[-which(report$Name %in% tree$tip.label),"Name"])
  
  
  report<-unique(report)
  
  colnames(report)[17:18]<-c("COVERAGE","IDENTITY")
  
  report$COVERAGE<-as.character(report$COVERAGE)
  report$IDENTITY<-as.character(report$IDENTITY)
  
  #ref_name<-sub(".*\\/", "", ref_name)
  
  report<-report[which(report$Name!="Name"),]
  genes_list<-grep("Gene",report$Gene.symbol,value = T,invert = T)
  
  report<-report[which(as.character(report$Gene.symbol) %in% genes_list),]
  
  report<-report[grep("fna|fasta",report$Name),]
  
  report$Name<-gsub(".fna","",report$Name)
  report$Name<-gsub(".fasta","",report$Name)
  
  # refname<-unique(report[-which(report$Name %in% tree$tip.label),"Name"])
  # 
  # report$Name<-ifelse(report$Name==refname,report$Contig.id,report$Name)
  # 
  report_res<-report[which(report$COVERAGE=="100.00" & report$IDENTITY=="100.00"),]
  report_res<-report_res[which(report_res$Element.type=="AMR" | report_res$Element.type=="STRESS"),]
  report_vir<-report[which(as.numeric(report$COVERAGE)>80 & as.numeric(report$IDENTITY)>80),]
  report_vir<-report_vir[which(report_vir$Element.type=="VIRULENCE"),]
  
  
  
  #report_summary<-rbind(report_res,report_vir)
  #pp<-table(report_summary$Name,report_summary$Gene.symbol)
  
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
  
  
  if (dim(vir_mtab1)[1]==0){
    summary_resvir<-res_mtab1
    summary_resvir[summary_resvir=="0"]<-NA
    
  } else if (dim(res_mtab1)[1]==0){
    summary_resvir<-vir_mtab1
    summary_resvir[summary_resvir=="0"]<-NA
    
  }else{
    
  
    vir_mtab1<-vir_mtab1[rownames(res_mtab1),]
    summary_resvir<-cbind(vir_mtab1,res_mtab1)
    summary_resvir[summary_resvir=="0"]<-NA
    
  }
  
  
  
  pp <- p %<+% rr + geom_tiplab(aes(fontface=length(rr$STRAIN_ID)*0.005),alpha=1,size=labsize*0.75,offset = 0.1,show.legend = FALSE)+
    theme_tree2()+theme(legend.text=element_text(size=length(rr$STRAIN_ID)*0.5),legend.title = element_text(size=length(rr$STRAIN_ID)*0.5))
  
  nclust<-length(rr$STRAIN_ID)*0.1
  
  rownames(rr)<-rr$STRAIN_ID
  rr$STRAIN_ID<-NULL
  
  suppressMessages(
    
    res_hm <- gheatmap(pp,summary_resvir, offset = 5, width=ncol(res_mtab1)*0.3/16, font.size=nclust, colnames_position= "top",
                       colnames_angle = 90, colnames_offset_y = -0.2, hjust = 0) +
      scale_y_continuous(limits=c(-1, NA))+
      scale_fill_manual(breaks=c("AMR/STRESS", "VIR"),
                        values=c("#6495ED", "firebrick"), name="Category",na.translate = F)
    
  )
  # rownames(rr)<-rr$STRAIN_ID
  # rr$STRAIN_ID<-NULL
  # 
  cluster_hm_levels <-res_hm + new_scale_fill()
  
  
  suppressMessages(
    
    
    hm_cluster <- gheatmap(cluster_hm_levels,rr, offset = 3, width=0.02,font.size = 5)+
      scale_fill_manual(values = color_vec,na.translate = F,name="Cluster")+
      scale_y_continuous(limits=c(-1, NA))+theme(legend.key.size = unit(nclust*0.5, 'cm'),
                                                 legend.spacing.y = unit(nclust*0.5,"cm"))
    
    
  )
  
  ggsave(hm_cluster, filename = "annotated_tree.svg",width = length(tree$tip.label)*1.5, height=length(tree$tip.label)*1.7 ,
         units = "cm",scale=1.2,limitsize = F)
  
  
  ###########  CLUSTER AND  NO RESISTANCE
  

  
} else if (file.exists("summary_resistance_virulence.txt")==FALSE) {
  
  
  pp <- p %<+% rr + geom_tiplab(aes(fontface=length(rr$STRAIN_ID)*0.005),alpha=1,size=labsize*0.75,offset = 0.1,show.legend = FALSE)+
    theme_tree2()+theme(legend.text=element_text(size=length(rr$STRAIN_ID)*0.5),legend.title = element_text(size=length(rr$STRAIN_ID)*0.5))
  
  nclust<-length(rr$STRAIN_ID)*0.1
  
  rownames(rr)<-rr$STRAIN_ID
  rr$STRAIN_ID<-NULL
  
  suppressMessages(
    
    
    hm_cluster <- gheatmap(pp,rr, offset = 2, width=0.02,font.size = labsize)+
      scale_fill_manual(values = color_vec,na.translate = F,name="Cluster")+
      scale_y_continuous(limits=c(-1, NA))+theme(legend.key.size = unit(nclust*0.5, 'cm'),
                                                 legend.spacing.y = unit(nclust*0.5,"cm"))
    
  )
  
  ggsave(hm_cluster, filename = "annotated_tree.svg",width = length(tree$tip.label)*1.5, height=length(tree$tip.label)*1.7 ,
         units = "cm",scale=1.2,limitsize = F)
  
  
  ######### NO CLUSTER / YES RESISTANCE
} else if (file.exists("summary_resistance_virulence.txt")) {
  
  
  report<-read.delim("summary_resistance_virulence.txt",sep='\t')
  
  refname<-unique(report[-which(report$Name %in% tree$tip.label),"Name"])
  
  
  report<-unique(report)
  
  colnames(report)[17:18]<-c("COVERAGE","IDENTITY")
  
  report$COVERAGE<-as.character(report$COVERAGE)
  report$IDENTITY<-as.character(report$IDENTITY)
  
  #ref_name<-sub(".*\\/", "", ref_name)
  
  report<-report[which(report$Name!="Name"),]
  genes_list<-grep("Gene",report$Gene.symbol,value = T,invert = T)
  
  report<-report[which(as.character(report$Gene.symbol) %in% genes_list),]
  
  report<-report[grep("fna|fasta",report$Name),]
  
  report$Name<-gsub(".fna","",report$Name)
  report$Name<-gsub(".fasta","",report$Name)
  
  # refname<-unique(report[-which(report$Name %in% tree$tip.label),"Name"])
  # 
  # report$Name<-ifelse(report$Name==refname,report$Contig.id,report$Name)
  # 
  report_res<-report[which(report$COVERAGE=="100.00" & report$IDENTITY=="100.00"),]
  report_res<-report_res[which(report_res$Element.type=="AMR" | report_res$Element.type=="STRESS"),]
  report_vir<-report[which(as.numeric(report$COVERAGE)>80 & as.numeric(report$IDENTITY)>80),]
  report_vir<-report_vir[which(report_vir$Element.type=="VIRULENCE"),]
  
  
  
  #report_summary<-rbind(report_res,report_vir)
  #pp<-table(report_summary$Name,report_summary$Gene.symbol)
  
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
  
  
  if (dim(vir_mtab1)[1]==0){
    summary_resvir<-res_mtab1
    summary_resvir[summary_resvir=="0"]<-NA
    
  } else if (dim(res_mtab1)[1]==0){
    summary_resvir<-vir_mtab1
    summary_resvir[summary_resvir=="0"]<-NA
    
  }else{
    vir_mtab1<-vir_mtab1[rownames(res_mtab1),]
    summary_resvir<-cbind(vir_mtab1,res_mtab1)
    summary_resvir[summary_resvir=="0"]<-NA
    
  }
  
  
  
  pp <- p %<+% rr + geom_tiplab(aes(fontface=length(rr$STRAIN_ID)*0.005),alpha=1,size=labsize*0.75,offset = 0.1,show.legend = FALSE)+
    theme_tree2()+theme(legend.text=element_text(size=length(rr$STRAIN_ID)*0.5),legend.title = element_text(size=length(rr$STRAIN_ID)*0.5))
  
  nclust<-length(rr$STRAIN_ID)*0.1
  
  rownames(rr)<-rr$STRAIN_ID
  rr$STRAIN_ID<-NULL
  
  suppressMessages(
    
    res_hm <- gheatmap(pp,summary_resvir, offset = 5, width=ncol(res_mtab1)*0.3/16, font.size=nclust, colnames_position= "top",
                       colnames_angle = 90, colnames_offset_y = -0.2, hjust = 0) +
      scale_y_continuous(limits=c(-1, NA))+
      scale_fill_manual(breaks=c("AMR/STRESS", "VIR"),
                        values=c("#6495ED", "firebrick"), name="Category",na.translate = F)
    
  )
  
  ggsave(res_hm, filename = "annotated_tree.svg",width = length(tree$tip.label)*1.5, height=length(tree$tip.label)*1.7 ,
         units = "cm",scale=1.2,limitsize = F)
  
  
  
  
  
}
