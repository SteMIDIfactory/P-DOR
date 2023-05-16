rm(list=ls())

suppressPackageStartupMessages({
  
  library(ggnewscale)
  library(ggplot2)
  library(svglite)
  library(dplyr)
  library(RColorBrewer)
  
})


args=commandArgs(trailingOnly = TRUE)


hospitalization_file <- args[1]
#cluster_file <- args[2]

hospitalization<-read.delim(hospitalization_file,sep='\t')

#hospitalization<-read.delim("metadata.txt",sep='\t')
###########################################################################
hospitalization$STRAIN_ID<-gsub(".fna","",hospitalization$STRAIN_ID)
hospitalization$STRAIN_ID<-gsub(".fasta","",hospitalization$STRAIN_ID)
hospitalization$PATIENT_ID<-as.character(hospitalization$PATIENT_ID)
hospitalization$WARD<-as.character(hospitalization$WARD)
hospitalization$WARD_DATE_ENTRY<-as.POSIXct(as.character(hospitalization$WARD_DATE_ENTRY))
hospitalization$WARD_DATE_EXIT<-as.POSIXct(as.character(hospitalization$WARD_DATE_EXIT))
###################################################################################
hospitalization$ISOLATION_DATE<-as.POSIXct(as.character(hospitalization$ISOLATION_DATE))
######################## MOST OCCURRENT BACTERIA
hospitalization$MATERIAL<-as.character(hospitalization$MATERIAL)
hospitalization$status<-ifelse(hospitalization$MATERIAL != "rectal_swab", "Infection", "Colonization")
hospitalization$STRAIN_ID<-as.character(hospitalization$STRAIN_ID)
############################################### MOST OCCURRENT WARD
hospitalization$WARD<-as.character(hospitalization$WARD)
most_ward<-head(sort(table(hospitalization$WARD),decreasing = TRUE), 5) #massimo 5
names_ward<-row.names(as.matrix(most_ward))
WARD2<-ifelse(hospitalization$WARD %in% names_ward, hospitalization$WARD, "Other_ward")
################################################ 

hospitalization<-cbind.data.frame(hospitalization, WARD2)

#############################################################################

cluster_file<-read.csv(list.files(pattern = "csv"),sep='')
cluster_file$cluster<-as.character(cluster_file$cluster)
cluster_file$STRAIN_ID<-as.character(cluster_file$STRAIN_ID)

  
out<-grep("BD",cluster_file$STRAIN_ID,value = T)

if (length(out)>0){
  
cluster_file<-cluster_file[-which(cluster_file$STRAIN_ID %in% out),]

freq_clust<-as.data.frame(table(cluster_file$cluster))
colnames(freq_clust)<-c("cluster","freq")

hosp_clust<-freq_clust[which(freq_clust$freq>1),"cluster"]
cluster_file<-cluster_file[which(cluster_file$cluster %in% hosp_clust),]

} else { 

cluster_file<-read.csv(list.files(pattern = "csv"),sep='')
cluster_file$cluster<-as.character(cluster_file$cluster)
cluster_file$STRAIN_ID<-as.character(cluster_file$STRAIN_ID)
}

########################################################################
mm<-merge(hospitalization,cluster_file,by = "STRAIN_ID",all = T)
clust<-select(mm,ISOLATION_DATE,PATIENT_ID, cluster)
clust<-clust[which(!is.na(clust$cluster)),]
moving<-unique(select(mm,PATIENT_ID,WARD_DATE_ENTRY, WARD_DATE_EXIT,WARD2))
moving<-moving[which(!is.na(moving$PATIENT_ID)),]
isol<-select(mm,ISOLATION_DATE,PATIENT_ID, status)
isol<-isol[which(!is.na(isol$PATIENT_ID)),]
###


ncols<-nrow(table(moving$WARD2))
pal_ward <- colorRampPalette(brewer.pal(5, "Set2"))(ncols)
ward<-pal_ward[1:length(levels(factor(moving$WARD2)))]
names(ward) <- levels(factor(moving$WARD2))


#pal_ward<-c("#e8177d","#7c44ad","#0072B2","#336a8f" ,"#2b9900" ,"#B3B3B3")
#ward<-pal_ward[1:length(levels(factor(moving$WARD2)))]
#names(ward) <- levels(factor(moving$WARD2))
#moving<-merge(moving,mm)

suppressWarnings(contact_network_plot<-ggplot(moving, start=WARD_DATE_ENTRY, end=WARD_DATE_EXIT, aes(WARD_DATE_ENTRY, PATIENT_ID)) +
                   theme_classic()+
                   geom_segment(data = moving,mapping = aes(color=WARD2, group=PATIENT_ID,
                   y=PATIENT_ID,xend=WARD_DATE_EXIT,yend=PATIENT_ID),alpha = .7, size=length(moving$PATIENT_ID)*0.7)+
                   scale_x_datetime(date_breaks='7 day',date_labels = "%Y-%m-%d")+
                   geom_point(data=isol, aes(ISOLATION_DATE,PATIENT_ID, shape=status),
                              #position=position_jitter(h=0.02),
                              alpha=1,size=length(moving$PATIENT_ID)*0.4,stroke = 0.5)+scale_shape_discrete(solid=F)+
                   labs(x = NULL, y = NULL, color = "ward")+scale_color_manual(values = ward)+
                   theme(axis.text.x=element_text(angle = 45, hjust = 1,size=length(moving$WARD_DATE_ENTRY)*0.9),
                         axis.text.y=element_text(hjust = 1,size=length(moving$PATIENT_ID)*0.9),
                         legend.text=element_text(size=length(moving$PATIENT_ID)*0.8),
                         legend.title =element_text(size=length(moving$PATIENT_ID)*1.2),
                         axis.title.x = element_text(hjust = 1,size=length(moving$PATIENT_ID)*0.7))+
                   labs(shape="condition")+ scale_linetype(guide = "none") +ggnewscale::new_scale_color() +
                   geom_line(data=clust, aes(ISOLATION_DATE,PATIENT_ID, group=cluster, linetype=cluster,color=cluster),
                             size=length(moving$PATIENT_ID)*0.08,alpha=1)+
                   scale_color_manual(values = brewer.pal(length(levels(factor(clust$cluster))), "Set2"))
                   
)



#adj<-ggplot_build(contact_network_plot)
var_height<-length(ggplot_build(contact_network_plot)$layout$panel_params[[1]]$y$get_labels())
var_width<-length(ggplot_build(contact_network_plot)$layout$panel_params[[1]]$x$get_labels())


ggsave(plot = contact_network_plot, filename = "contact_network_plot.svg",
       width = var_width*3, height=var_height*3 , units = "cm",limitsize = F,scale = 1)
