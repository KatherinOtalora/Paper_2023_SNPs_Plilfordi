############### R PIPELINE SNPS RIQUEZA ALELICA ############

### Abrir otra ventana en R #######
#Ir a la terminal y ahi ir hasta la carpeta de aplicaciones de Mac luego ejecutar

open -n /Applications/RStudio.app    # Para abir varias ventanas de R, ejecutar esto desde la terminal de tu MAC.

#Instalar paquetes necesarios para el analisis de SNPs

install.packages("vcfR")
install.packages("poppr")
install.packages("ape")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("pillar")
install.packages("adegenet")
install.packages("dartR") #No se pudo instalar en R studio
install.packages("ggplot2")
install.packages("pegas")
install.packages("hierfstat") #No se pudo instalar en R studio pero si en R version 4.1
install.packages("devtools")
install_github("green-striped-gecko/dartR") #No se pudo instalar en Rstudio

if (!requireNamespace("BiocManager", quietly = TRUE))
  +     install.packages("BiocManager")
browseVignettes("SNPRelate")
install.packages("ggplot2")
install.packages("factoextra")

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")
install.packages("tidyr")

#Lllamar las librerias necesarias para el analisis de SNPs

library(vcfR)
library(poppr)
library(ape)
library(pheatmap)
library(RColorBrewer)
library(pillar)
library(adegenet)
library(dartR)
library(adegenet)
library(hierfstat)
library(devtools)
library(pegas)
library(ggplot2)
library(factoextra)
library(tidyr)

############# OJO AQUI EL ARHIVO QUE SE USA ES EL QUE NO TIENE EL LD

#Nos ubicamos en la carpeta local desde donde vamos a trabajar, hay que asegúrese de que está en la carpeta correcta con los archivos descargados disponibles.

setwd("/Users/katherin_otalora/Documents/R_analysis/GATK_2022/allSamples/allSamples_0.90")

#El output de LD from PLINK tiene los individuos con doble nombre y hay que dejarlos con su nombre individual, se puede cambiar con texwangler esto es por el --double-id - told plink to duplicate the id of our samples

# Asignamos  nuestros datos .vcf filtrados por LD a una variable nueva con nombre podarcislilfordi.VCF

podarcislilfordi.VCF <- read.vcfR("allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.vcf")

podarcislilfordi.VCF

#Los datos VCF no suelen incluir ningún tipo de información sobre la población. Tenemos que cargar estos datos por separado del archivo delimitado por texto que descargamos anteriormente y que incluye el ID de las muestras y el estado del que se obtuvieron las muestras. Carguemos el archivo population_data.gbs.txt en la memoria utilizando la función read.table():

pop.data <- read.table("popfile_3.txt", sep = "\t", header = TRUE)

#Ahora podemos comprobar que todas las muestras del VCF y del archivo de datos de la población están incluidas correctamente:

#Compruebe los nombres en los dos archivos
colnames(podarcislilfordi.VCF@gt)
colnames(podarcislilfordi.VCF@gt)[-1]#the option [-1] eliminates the "Format" cell
pop.data$AccessID

#Por último, compruebe la identidad
#Tiene que dar TRUE para corroborar que mi data set VCF y el pop.data son iguales si da FALSE hay que revisar si el pop.data esta mal

all(colnames(podarcislilfordi.VCF@gt)[-1] == pop.data$AccessID)

#Si no coiciden hay que ir a ver los individuos del vcf y corroborar que sean los mismos de popfile, aveces son diferentes porque hemos eliminado individuos con filtros anteriores

#Creamos un objeto genlight con vcfR2genlight
#El paquete vcfR contiene múltiples funciones para convertir los datos en otros formatos. Para nuestro propósito particular queremos convertir el objeto vcfR en un objeto genlight. Para ello podemos utilizar la función vcfR2genlight:

gl.podarcislilfordi <- vcfR2genlight(podarcislilfordi.VCF)

#Al transformar el objeto se muestra una advertencia que nos indica que hay algunos loci con más de dos alelos. Muchas de las funciones que utilizaremos han sido escritas bajo el supuesto de un modelo bialélico. Este modelo bialélico restringe todos los loci a sólo dos alelos para simplificar algunos cálculos. El objeto genlight sólo admite loci con no más de dos alelos. La función vcfR2genlight organiza los datos para filtrar los loci que no son bialélicos, devolviendo un file que sólo contiene loci con dos alelos. La advertencia es para asegurarnos de que somos conscientes de que esta acción ha tenido lugar.

#Además, se nos pide que especifiquemos la ploidía del organismo para poder calcular algunas métricas genéticas de la población. P.lilfordi es un organismo diploide, por lo que especificaremos una ploidía de dos. Todos los objetos genlight tienen opciones de ploidía, en las que el usuario puede especificar la ploidía de cada muestra individual, o una vez para toda la población. Podemos asumir que cada muestra de P.lilfordi es diploide y especificaremos una ploidía de 2 para toda la población. hay que tener  en cuenta que aunque un objeto genlight puede admitir individuos de diferente ploidía, dentro de cada individuo todos los loci deben ser de la misma ploidía.

ploidy(gl.podarcislilfordi) <- 2

#Nuestra pregunta biológica para P.lilfordi requiere poblaciones predeterminadas. Podemos añadirlas al objeto genlight como parte de la opcion pop (población). Para especificar la población, añadimos la columna Locality de nuestro set de datos pop.data a la opcion pop de nuestro objeto genlight:

pop(gl.podarcislilfordi) <- pop.data$Locality

#Ahora tenemos un objeto genlight de datos VCF filtrados:

gl.podarcislilfordi

#TRANSFORMAR GENLIGHT A GENIND OBJECT esto se realiza ya que varios analisis necesitan este input para poder ejecutarse

gnd.podarcislilfordi <-gl2gi(gl.podarcislilfordi, probar = FALSE, verbose = NULL) #Para convertir genlight en genind ... Converts a genlight object to genind object
save(gnd.podarcislilfordi, file = "gnd.podarcislilfordi.RData")

gnd.podarcislilfordi

############################# RIQUEZA ALELICA ###########################
#Riqueza alélica (AR): Número promedio de alelos por locus independiente del tamaño de muestra.

#Riqueza alelica por locus
Allelic_stats_richness <- allelic.richness(genind2hierfstat(gnd.podarcislilfordi),min.n=8,diploid=TRUE)$Ar
#Allelic_stats_richness <- allelic.richness(genind2hierfstat(gnd.podarcislilfordi),min.n=NULL,diploid=TRUE)$Ar
#Allelic_stats_richness <- allelic.richness(genind2hierfstat(gnd.podarcislilfordi))$Ar
Allelic_stats_richness_datos_limpios <- na.omit (Allelic_stats_richness) #Para eliminar los NA de las estadísticas
write.table(Allelic_stats_richness_datos_limpios , "Allelic_stats_richness_minn8_per_locus_allsamples.txt", quote=FALSE, row.names=TRUE)

#Riqueza alelica por poblacion

media_addaiagran <- round (mean(Allelic_stats_richness_datos_limpios$ADGR) %>% round(digits = 7)) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$ADGR [is.na(Allelic_stats_richness$ADGR)] <- media_addaiagran #Remplazar los NA por su media en la poblacion

media_aire <- round (mean(Allelic_stats_richness_datos_limpios$AIRE) %>% round(digits = 7)) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$AIRE [is.na(Allelic_stats_richness$AIRE)] <- media_aire #Remplazar los NA por su media en la poblacion

media_cabrerah <- round (mean(Allelic_stats_richness_datos_limpios$CABRERAH, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$CABRERAH [is.na(Allelic_stats_richness$CABRERAH)] <- media_cabrerah #Remplazar los NA por su media en la poblacion

media_cabrerai <- round (mean(Allelic_stats_richness_datos_limpios$CABRERAI, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$CABRERAI [is.na(Allelic_stats_richness$CABRERAI)] <- media_cabrerai #Remplazar los NA por su media en la poblacion

media_colom <- round (mean(Allelic_stats_richness_datos_limpios$COLOM, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$COLOM [is.na(Allelic_stats_richness$COLOM)] <- media_colom #Remplazar los NA por su media en la poblacion

media_colomer <- round (mean(Allelic_stats_richness_datos_limpios$COLOMER, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$COLOMER [is.na(Allelic_stats_richness$COLOMER)] <- media_colomer #Remplazar los NA por su media en la poblacion

media_dragonera <- round (mean(Allelic_stats_richness_datos_limpios$DRAGONERA, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$DRAGONERA [is.na(Allelic_stats_richness$DRAGONERA)] <- media_dragonera #Remplazar los NA por su media en la poblacion

media_esclts <- round (mean(Allelic_stats_richness_datos_limpios$ESCLTS, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$ESCLTS [is.na(Allelic_stats_richness$ESCLTS)] <- media_esclts #Remplazar los NA por su media en la poblacion

media_escurt <- round (mean(Allelic_stats_richness_datos_limpios$ESCURT, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$ESCURT [is.na(Allelic_stats_richness$ESCURT)] <- media_escurt #Remplazar los NA por su media en la poblacion

media_foradada <- round (mean(Allelic_stats_richness_datos_limpios$FORADADA, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$FORADADA [is.na(Allelic_stats_richness$FORADADA)] <- media_foradada #Remplazar los NA por su media en la poblacion

head (Allelic_stats_richness) #Para ir viendo los cambios que se van haciendo

media_guardia <- round (mean(Allelic_stats_richness_datos_limpios$GUARDIA, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$GUARDIA [is.na(Allelic_stats_richness$GUARDIA)] <- media_guardia #Remplazar los NA por su media en la poblacion

media_moltona <- round (mean(Allelic_stats_richness_datos_limpios$MOLTONA, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$MOLTONA [is.na(Allelic_stats_richness$MOLTONA)] <- media_moltona #Remplazar los NA por su media en la poblacion

media_porros <- round (mean(Allelic_stats_richness_datos_limpios$PORROS, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$PORROS [is.na(Allelic_stats_richness$PORROS)] <- media_porros#Remplazar los NA por su media en la poblacion

media_porroscaballeria <- round (mean(Allelic_stats_richness_datos_limpios$PRCV, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$PRCV [is.na(Allelic_stats_richness$PRCV)] <- media_porroscaballeria #Remplazar los NA por su media en la poblacion

media_rovells <- round (mean(Allelic_stats_richness_datos_limpios$RAVL, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$RAVL [is.na(Allelic_stats_richness$RAVL)] <- media_rovells #Remplazar los NA por su media en la poblacion

media_rei <- round (mean(Allelic_stats_richness_datos_limpios$REI, digits = 7)) %>% round(digits = 7) #Para hacer un variable que calcule la media por poblacion con los datos limpios
Allelic_stats_richness$REI [is.na(Allelic_stats_richness$REI)] <- media_rei #Remplazar los NA por su media en la poblacion

Per_population <- Allelic_stats_richness %>% apply(MARGIN = 2, FUN = mean) %>% round(digits = 6) # Para ver la riqueza alelica por poblaciones
Per_population
write.table(Per_population , "Allelic_stats_richness_minn8_per_population_allsamples.txt", quote=FALSE, row.names=TRUE)
#min.n=8 Se pone 8 por que si quiero 8 alelos como mínimo, necesitaré por lo menos 4 individuos (2 alelos x 4 indiv=8 allelos)
#min.n es el numero de allelos y min.all en el output es el numero de allelos usado para la rarefacción.
#Si tu le pones min.n=8, son 8 alelos minimo para la rarefacción, lo que quiere decir, que siendo los individuos diploides (cada indiv tiene como max 2 allelos por locus) que, si quiero 8 alelos como mínimo, necesitaré por lo menos 4 individuos (2 alelos x 4 indiv=8 allelos). 
#Si pongo min.n= 16, necesitaré 8 individuos como minimo, entonces para Dragonera y Colomer no puedo calcular la AR ya que solo tienen 4 indiv

#ADGR      AIRE  CABRERAH  CABRERAI     COLOM   COLOMER DRAGONERA    ESCLTS    ESCURT 
#1.349890  1.395833  1.430697  1.410910  1.394326  1.357958  1.437806  1.317870  1.348026 
#FORADADA   GUARDIA   MOLTONA    PORROS      PRCV      RAVL       REI 
#1.337569  1.422111  1.413981  1.244472  1.344511  1.342078  1.373006 
