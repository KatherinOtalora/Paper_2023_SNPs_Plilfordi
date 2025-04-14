################### Identifying candidate genes with GFF file EXTSAMPLES #######################

#Ahora que tenemos SNPs que creemos que pueden estar bajo selección, podemos dar el siguiente paso para identificar cuáles son los genes que se encuentran cerca de esos SNPs. Haremos esto completamente en R.
#Entrar a la unicacion de la carpeta:

setwd ("/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Filtro_LOG_PO_FST_extsamples")

#En primer lugar vamos a leer nuestros outliers comunes entre Bayescan y PCAdapt
#Cargamos las librerias de nuestro interes

library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

#Abra el archivo de salida de los outliers comunes entre pcadapt y bayescan con la extensión "st_statisc.txt"**.

setwd ("/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Filtro_LOG_PO_FST_extsamples")

#Importante los nombres de los SNP hay que serparlos por un lado los cromosomas y por otro la posicion rPodLil1.1_Chr1_46895223 poner un TAB en el _ en el .txt y no poner nombres de las columnas (lo haremos mas adelante)
#Leemos nuestro archivo ya editado

bayescan_and_pcadapt= read.table("/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Filtro_LOG_PO_FST_extsamples/outliers_comunes_extSamples_1_100_K=10_fst_statisc.txt")
head (bayescan_and_pcadapt,10)

#La primera columna del marco de datos de bayescan es el ID del SNP (SNP y POS en nuestro caso).
#Las tres columnas siguientes (prob, log10(P0), y qval) están relacionadas con la prueba de adaptación local considerando el logaritmo de las probabilidades posteriores - log10(PO) - y el valor q para el modelo con selección.
#La quinta columna da el tamaño del efecto específico del locus (parámetro alfa). La última proporciona el FST específico del locus promediado en todas las poblaciones.

#Cambiar el nombre de las columnas.

colnames(bayescan_and_pcadapt)=c("CHR","POS","PROB","LOG_PO","Q_VALUE","ALPHA","FST")
head (bayescan_and_pcadapt,10)

#Ahora vamos a filtrar por FST > 0.80 Y LOG_PO==1000

OUTLIERS_filtred_LOGPO <- bayescan_and_pcadapt %>% arrange(LOG_PO) %>% filter(LOG_PO==1000)
head (OUTLIERS_filtred_LOGPO,10)

OUTLIERS_filtred_LOGPO_FST <- OUTLIERS_filtred_LOGPO %>% arrange(FST) %>% filter(FST > 0.80)
head (OUTLIERS_filtred_LOGPO_FST,10)

dim (OUTLIERS_filtred_LOGPO_FST)

write.table(OUTLIERS_filtred_LOGPO_FST, "outliers_filtred_LOGPO_FST_extsamples_withoutLD_1_100_K=10_.txt", sep="\t", quote=FALSE)
#Sacamos los nombres de los outliers para usar el gff y buscar a que proteina codifican
###write.table(new_gff$chrnew_gff$pos, "EXP.txt", sep="\t", quote=FALSE)


#IR A VCFTOOLS Y EXTRAER LOS OUTLIERS filtrados por LOGPO_FST
# Hacer un archivo con los outrliers de bayescan con la informacion del cromosoma un TAB y luego la posicion
#En vcftools la opcion positions me permite extraer los snps que me interesan

cd /Users/katherin_otalora/Documents/GFF_Genome_Plilfordi//Filtro_LOG_PO_FST_extsamples
vcftools --vcf extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf --positions outliers_filtred_LOGPO_FST_extamples_withoutLD_1_100_K=10_positions.txt --recode --recode-INFO-all --out outliers_filtred_LOGPO_FST_comunes_bayescan_pcadapt_extsamples_withoutLD_1_100_K=10

grep -v # outliers_filtred_LOGPO_FST_comunes_bayescan_pcadapt_extsamples_withoutLD_1_100_K=10.recode.vcf | wc -l
#190

#FIN
