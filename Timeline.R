rm(list=ls())
library(plyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(reshape2)
library(wesanderson)
#input

args=commandArgs(trailingOnly = TRUE)


ricoveri_file <- args[1]


ricoveri<-read.delim(ricoveri_file,sep='\t')
colnames(ricoveri)[1]<-"ID_registry"
colnames(ricoveri)[5]<-"ward_name"
colnames(ricoveri)[3]<-"ward_date_entry"
colnames(ricoveri)[4]<-"ward_date_exit"
###########################################################################
ricoveri$ID_registry<-as.character(ricoveri$ID_registry)
ricoveri$ward_name<-as.character(ricoveri$ward_name)
ricoveri$ward_date_entry<-as.POSIXct(as.character(ricoveri$ward_date_entry))
ricoveri$ward_date_exit<-as.POSIXct(as.character(ricoveri$ward_date_exit))
###################################################################################
ricoveri$ISOLATION_DATE<-as.POSIXct(as.character(ricoveri$ISOLATION_DATE))
######################## MOST WANTED BACTERIA
ricoveri$MATERIAL<-as.character(ricoveri$MATERIAL)
most_material<-head(sort(table(ricoveri$MATERIAL),decreasing = TRUE), 5) #massimo 5
names_material<-row.names(as.matrix(most_material))
sample_material2<-ifelse(ricoveri$MATERIAL %in% names_material, ricoveri$MATERIAL, "Other_material")
### OR
ricoveri$status<-ifelse(ricoveri$MATERIAL != "rectal_swab", "Infection", "Colonization")

############################################### MOST WANTED WARD
ricoveri$ward_name<-as.character(ricoveri$ward_name)
most_ward<-head(sort(table(ricoveri$ward_name),decreasing = TRUE), 5) #massimo 5
names_ward<-row.names(as.matrix(most_ward))
ward_name2<-ifelse(ricoveri$ward_name %in% names_ward, ricoveri$ward_name, "Other_ward")
################################################ aggiungo colonne col vettore rinominato max 5 + other al dataframe
ricoveri<-cbind.data.frame(ricoveri, ward_name2)
#setto una palette specifica per i colori in modo da aver sempre gli stessi abbinati
###################### aggiungo colori reparti
mostwward<-sort(table(ricoveri$ward_name2)) #prendo i primi 5 + other definiti dal dataframe precedente
names_mostwward<-row.names(as.matrix(mostwward)) #nomi dei reparti + other
pal_ward<-c("#336a8f", "#e8177d", "#8E44AD","#2b9900","#0072B2","#484848")

dd1<-cbind.data.frame(names_mostwward,pal_ward)
#mettere reparti come rownames e richiamare il dataframe per associare il colore in tutte le value
colnames(dd1)<-c("ward_name2", "ward_colors")
#row.names(dd1)<-dd1$ward_name2
ward_color_vector <- dd1[as.matrix(ricoveri$ward_name2), "ward_colors"]
ricoveri<-cbind.data.frame(ricoveri, ward_color_vector)
ward_color_vector<-as.character(ward_color_vector)
ricoveri$ward_color_vector<-ward_color_vector
ward_lvls  <- levels(ricoveri$ward_name2)
ward_cols  <- setNames(hcl(h=seq(length=length(ward_lvls))),ward_lvls)
ward_cols[ricoveri$ward_name2]<-ward_color_vector

#############################################################################

timelineplot <- ggplot(ricoveri, start=ward_date_entry, end=ward_date_exit, aes(ward_date_entry, ID_registry)) +
  theme_light()+
  geom_segment(aes(color=ward_name2, group=ID_registry,
                   y=ID_registry,xend=ward_date_exit,yend=ID_registry),alpha = 0.75, size=6) +
  scale_x_datetime(date_breaks='7 day',date_labels = "%Y-%m-%d")+theme(axis.text.x=element_text(angle=90,hjust = 1, size = 7))+
  geom_point(data=ricoveri, aes(ISOLATION_DATE,ID_registry, shape=status),
             position=position_jitter(h=0.02),
             alpha=1,size=3,stroke = 0.5)+labs(shape="Condition")+scale_shape_discrete(solid=F)+
  labs(x = 'Hospitalization', y = NULL, color = "Ward")+
  theme(axis.text.x=element_text(angle = 45, hjust = 1,size=10),
        axis.text.y=element_text(hjust = 1,size=12),
        legend.text=element_text(size=12))


#
#
#var_height<-length(unique(ricoveri$ID_registry))*0.03
#

ggsave(timelineplot, filename = "timeline.png",width = 40, height =40 , units = "cm")
# ##############################################################################################################
