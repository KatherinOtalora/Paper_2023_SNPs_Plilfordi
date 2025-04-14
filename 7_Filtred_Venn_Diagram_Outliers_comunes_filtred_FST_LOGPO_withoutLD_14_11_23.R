####################################################### VENN DIAGRAM ###################################

#Cargar las librerias

install.packages("VennDiagram")
library(VennDiagram)
install.packages("ggVennDiagram")
library(ggVennDiagram)
install.packages("sf")
library(sf)

#Guardar cada set de outliers en un vector independiente

setwd("/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi")


outliers_comunes_filtred_intsample <- read.table ("/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/names_outliers_filtred_LOGPO_FST_intsamples_withoutLD_1_100_K=8_1.txt")
outliers_comunes_filtred_extsample <-read.table ("/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/names_outliers_filtred_LOGPO_FST_extamples_withoutLD_1_100_K=10_1.txt")
outliers_comunes_filtred_allsample_0.90 <- read.table ("/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/names_outliers_filtred_LOGPO_FST_allsamples_withoutLD_1_100_K=15_1.txt")

x_int <- unlist(outliers_comunes_filtred_intsample)
y_ext <- unlist(outliers_comunes_filtred_extsample)
z_all <- unlist(outliers_comunes_filtred_allsample_0.90)


length(x_int)
length(y_ext)
length(z_all)


#Crea un diagrama de VENN

--------------------------------------------------------------
library(ggVennDiagram)

x <- list(Ext=y_ext, Comb=z_all, Int=x_int);x
y <- list(Ext=y_ext,Int=x_int);x

ggVennDiagram (x[1:3]) # 3 data set
ggVennDiagram (x[1:2])  #Ext and Comb
ggVennDiagram (x[2:3]) #Comb and Int
ggVennDiagram (y[1:2]) #Ext and Int

#--------------------------------------------------------------
outliers_comunes_ext_all <- intersect(y_ext, z_all)
outliers_comunes_ext_all
write.table(outliers_comunes_ext_all, "outliers_comunes_ext_all_filtred_FST_LOPPo.txt", row.names=F, quote=F)

outliers_comunes_int_all <- intersect(x_int, z_all)
outliers_comunes_int_all
write.table(outliers_comunes_int_all, "outliers_comunes_int_all_filtred_FST_LOPPo.txt", row.names=F, quote=F)

outliers_comunes_int_ext <- intersect(x_int, y_ext)
outliers_comunes_int_ext
write.table(outliers_comunes_int_ext, "outliers_comunes_int_ext_filtred_FST_LOPPo.txt", row.names=F, quote=F)
