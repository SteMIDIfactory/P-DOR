library(ape)
library(writexl)
#install.packages("adegenet")

#setwd("/home/ghepard/progetto_CHIRURGIE_ERCP/filogenesi/")

args=commandArgs(trailingOnly = TRUE)
setwd(getwd())
snp_alignment <- args[1]


seq<-read.dna("SNP_alignment.core.fasta", format="fasta")
dist_matrix<-as.matrix(dist.dna(seq, "N"))

##########
outbreak_genomes<-as.vector(grep("DB_",rownames(dist_matrix),value = TRUE, invert = TRUE))
external_genomes<-as.vector(grep("DB_",rownames(dist_matrix),value = TRUE))
ref_genomes<-as.vector(grep("REF_",rownames(dist_matrix),value = TRUE))
outbreak_genomes<-outbreak_genomes[- which(outbreak_genomes %in% ref_genomes)]

snp_tab<-dist_matrix[outbreak_genomes,outbreak_genomes] #relativa ai genomi di outbreak di studio
out_snp_tab<-dist_matrix[external_genomes,external_genomes] #relativa ai genomi di outbreak di studio
##########

cairo_ps("all.eps", height = 12, width = 15, pointsize = 12, fallback_resolution = 300)
h0<-hist(dist_matrix,breaks=200, xlim = range(0, max(dist_matrix)))
#print(h0)
dev.off()

cairo_ps("patric.eps", height = 12, width = 15, pointsize = 12, fallback_resolution = 300)
h0<-hist(out_snp_tab,breaks=200, xlim = range(0, max(out_snp_tab)))
#print(h0)
dev.off()


cairo_ps("outbreak.eps", height = 12, width = 15, pointsize = 12, fallback_resolution = 300)
h0<-hist(snp_tab,breaks=200, xlim = range(0, max(snp_tab)))
#print(h0)
dev.off()


#Histogram Colored (blue and red)

#hist(h0$counts,col=rgb(1,1,0,0.5), main="Overlapping Histogram", xlab="SNPs distance",
 #    breaks=1000, xlim=range(0,100),ylim=range(0,100))

#hist(h1$counts, col=rgb(1,0,0,0.5), main="Overlapping Histogram", xlab="SNPs distance",
#     breaks=1000, xlim=range(0,100),ylim=range(0,100))

#hist(h2$counts, col=rgb(0,0,1,0.5), breaks=100, xlim=range(0,500),ylim=range(0,500),add=T)
#box()



write_xlsx(as.data.frame(dist_matrix),"all_snp_distance_tab.xlsx")
write_xlsx(as.data.frame(out_snp_tab),"patric_snp_distance_tab.xlsx")
write_xlsx(as.data.frame(snp_tab),"outbreak_snp_distance_tab.xlsx")

#dens<-density(out_snp_tab)
#plot(dens)
#dens<-density(snp_tab)
#plot(dens)
#write.xlsx(dist_matrix,"total_snp_distance_tab.xlsx")


