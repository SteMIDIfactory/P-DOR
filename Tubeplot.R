rm(list=ls())


library(RColorBrewer)
library(ggplot2)
library(ggstance)
library(reshape2)
library(data.table)
library(wesanderson)
library(scales)

###
args=commandArgs(trailingOnly = TRUE)


ricoveri_file <- args[1]
ricoveri<-read.delim(ricoveri_file,sep='\t')
colnames(ricoveri)[1]<-"ID"
#NUOVA FORMATTAZIONE
#############################################################
ricoveri$ID<-as.character(ricoveri$ID)
#ricoveri$INIZIORIC<-as.POSIXct(strptime(ricoveri$INIZIORIC,"%Y/%m/%d"))
Sys.setenv(TZ='Europe/Rome')
ricoveri$INIZIORIC<-as.POSIXct(as.character(ricoveri$INIZIORIC))
ricoveri$FINERIC<-as.POSIXct(as.character(ricoveri$FINERIC))
ricoveri$ISOLATION_DATE<-as.POSIXct(as.character(ricoveri$ISOLATION_DATE))
ricoveri$REPARTO<-as.character(ricoveri$REPARTO)
ricoveri$MATERIAL<-trimws(as.character(ricoveri$MATERIAL), "l")
ricoveri$MATERIAL<-trimws(as.character(ricoveri$MATERIAL),"r")
ricoveri$STATUS<-ifelse(ricoveri$MATERIAL == "rectal_swab", "Colonized","Infected")
###ifelse n paz = 1 usa geom line
one_ward<-ricoveri[which(table(ricoveri$ID)==1),]
###predo solo gli id dei pazienti con un singolo ingresso
ward<-as.data.frame(table(ricoveri$ID))
one_ward_id<-ward[which(ward$Freq==1),"Var1"]
'%nin%' <- Negate('%in%')
multiple_ward<-ricoveri[which(ricoveri$ID %nin% one_ward_id),]
one_ward<-ricoveri[which(ricoveri$ID %in% one_ward_id),]

### PRENDO L'ULTIMO INTERVALLO DI RICOVERO PER FARCI IL GEOM_SEGM E COMPLETARE IL PLOT
library(tidyverse)
end_ward<-multiple_ward %>%
  group_by(ID) %>%
  slice(tail(row_number(), 1))
end_ward<-as.data.frame(end_ward)

### INGANNO GEOM_STEP METTENDO UN ULTERIORE RICOVERO DUMMY PER FARGLI TIRARE LA LINEA FINALE
end_ward1<-end_ward
end_ward1$INIZIORIC<-end_ward1$FINERIC
end_ward1$STATUS<-NA
end_ward1$ISOLATION_DATE<-NA
end_ward1$MATERIAL<-NA

enddf<-rbind(multiple_ward,end_ward1)
#options(ggrepel.max.overlaps = Inf)
tube<-ggplot(enddf,aes(INIZIORIC,REPARTO, color=ID,group=ID))+theme_light()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1,size=10),
        axis.text.y=element_text(hjust = 1,size=12),
        legend.text=element_text(size=12),legend.title = element_blank())+
  geom_step(position = position_dodge2v(height = 0.35),direction = "hv",size=1,alpha=0.6)+
    geom_segment(data = one_ward,aes(x = INIZIORIC,xend=FINERIC,yend=REPARTO,group=ID),size=1,alpha=0.6,
                 position = position_jitter(height = 0.15))+
               #geom_text_repel(data = enddf,aes(ISOLATION_DATE,REPARTO,label=STATUS),box.padding = 2.5,direction= "y")+
               #geom_point(data = enddf,shape=23,size=2,position = position_dodge2v(height = 0.15, preserve = "total"))+
               scale_x_datetime(breaks=date_breaks("7 day"),date_labels = ("%Y-%m-%d"),timezone = "CEST")+
                  labs(x=NULL,y=NULL)


ggsave(tube, filename = "tube-map_plot.png",width = 10, height = 10 , units = "cm")
