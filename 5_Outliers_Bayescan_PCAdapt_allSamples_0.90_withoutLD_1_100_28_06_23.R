##################### Outliers WITHOUTLD ALLSAMPLES 0.90 K=15#############################

#BAYESCAN
#BayeScan http://cmpg.unibe.ch/software/BayeScan/ es un programa de línea de comandos cuyo objetivo es identificar loci candidatos putativos bajo selección natural a partir de datos genéticos, utilizando las diferencias en las frecuencias alélicas entre grupos especificados. Los grupos pueden establecerse utilizando los lugares de muestreo o las unidades genéticas cuando se investiga la estructura de la población.
#BayeScan se basa en el modelo Multinomial-Dirichlet.
#Este programa puede definir tres categorías de loci candidatos putativos:

#bajo selección diversificadora
#bajo selección equilibradora
#bajo neutralidad
#Para cada locus, BayeScan calcula una probabilidad posterior (Posterior odds) - disponible a través del parámetro pr_odds - para el modelo que incluye la selección. Estas probabilidades posteriores indican cuánto más probable es el modelo con selección en comparación con el modelo neutral. Por ejemplo, una pr_odds de 10 significa que hay una probabilidad de 1 en 10 de que un marcador esté bajo selección. Este número sería demasiado alto si se considera un conjunto de datos con hasta 10.000 #marcadores.

#En el contexto de las pruebas múltiples con un gran número de marcadores (hasta 10.000), ejecute BAYESCAN con los parámetros adecuados, como se recomienda en Whitlock y Lotterhos (2015) https://www.jstor.org/stable/10.1086/682949?seq=1.

#Para ello, debe tener en cuenta el número de loci de su conjunto de datos. También puede consultar el ejercicio de BayeScan https://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/bayescan-exercise/ para saber más sobre cómo interpretar los archivos y resultados de BayeScan.

#Desde la terminal:

cd /Users/katherin_otalora/Documents/PGDSpider_2.1.1.5

java -Xmx1024m -Xms512m -jar PGDSpider2.jar

#Se abrira una ventana emergente java interactiva para convertir los archivos, vamos a transformar un VCF to GENEPOP (Tenemos que organizarlos y demiliartlos con el POP antes de cada nueva poblacion)
#Luego transformamos los datos de GENEPOP (organizados y con el POP delimitando cada poblacion) y se pasan a GESTE BAYESCAN en PGDSpider (Se corrio en MAC "(base) katherin_otalora@KatherilorasMBP PGDSpider_2.1.1.5 % java -Xmx1024m -Xms512m -jar PGDSpider2.jar")
#Ponerle al final del archivo .txt para poder leerlo en BayeScan

############################ Run BayeScan in a terminal ###################################

#Entramos en la carpeta

brew install gcc   #Instalar gcc para ejecutar el Bayescan

#Para solucionar el error de "lib"
sudo cp /Users/katherin_otalora/Downloads/libgcc_s.1.dylib  /usr/local/lib/libgcc_s.1.dylib
sudo cp /Users/katherin_otalora/Downloads/libgomp.1.dylib  /usr/local/lib/libgomp.1.dylib

#Entramos a nuestra carpeta

cd /Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/allSamples_0.90_withoutLD_bayescan_1_100_K=15

#Corremos el Bayescan con 1:100
#Testeamos: 20 short pilots runs with 5000 integration, burn in set to 5 x 104 and thinning interval 10 and prior odds of 100.

/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/BayeScan2.1_macos64bits /Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/allSamples_0.90_withoutLD_bayescan_1_100_K=15/allSamples_0.90_filtredwithoutLD_bayescan_1_100.txt -od ./ -threads 4 \ -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100

#La primera línea consiste en las opciones básicas para asegurar que el programa se ejecuta correctamente - es decir, dónde está el archivo de entrada, en qué directorio escribir la salida (./ significa el directorio actual) y el número de hilos o núcleos que podemos utilizar (4 en este caso).

#La segunda línea especifica las opciones para los análisis y la cadena de Markov Monte Carlo y son las siguientes

#-n es el número de iteraciones que queremos muestrear del MCMC (5000 en este caso)
#-thin es el intervalo de thin, fijado en 10; esto significa que dejamos que el MCMC se ejecute durante 10 pasos entre cada muestra.
#-nbp es el número de ejecuciones piloto que Bayescan ejecuta para ajustar los parámetros del MCMC.
#-pilot es la duración de cada una de las ejecuciones piloto, 5000 iteraciones aquí.
#-burn es la longitud del burnin (es decir, la parte del MCMC que se descarta), aquí es de 50 000
#prior odds es un valor importante porque tiene un efecto sobre los falsos positivos: usamos 100 para ser conservadores, 10 es menos conservador

#Luego de correr saldran varios outputs en la misma carpeta donde se ubico.

###################### Utilice R para identificar los valores atípicos de los análisis de BayeScan ########################

install.packages("vcfR")
install.packages ("hierfstat")
install.packages ("adegenet")
install.packages ("ggplot2")
install.packages ("radiator")

library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

#Abra el archivo de salida de **BayeScan con la extensión "_fst.txt "**.

setwd("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/allSamples_0.90_withoutLD_bayescan_1_100_K=15")

bayescan=read.table("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/allSamples_0.90_withoutLD_bayescan_1_100_K=15/allSamples_0.90_filtredwithoutLD_bayescan_1_100_fst.txt")

#La primera columna del marco de datos de bayescan es el ID del SNP.
#Las tres columnas siguientes (prob, log10(P0), y qval) están relacionadas con la prueba de adaptación local considerando el logaritmo de las probabilidades posteriores - log10(PO) - y el valor q para el modelo con selección.
#La quinta columna da el tamaño del efecto específico del locus (parámetro alfa). La última proporciona el FST específico del locus promediado en todas las poblaciones.

#Descargue la lista de SNPs en el orden correcto. El formato .geste tiene los SNPs en el mismo orden que el vcf utilizado para producir el formato .geste. Por lo tanto, puede utilizar este comando en bash para extraer la tercera columna que contiene la información de identificación de cada SNP en su vcf.
#Importante para el archivo orginal de vcf hay que cortar la columna 1 y 2 y luego juntarlas porque en la columna tres no hay informacion solo un . (Se juntan en algun editor de texto con _ entre medio de la culmna 1 y 2)

grep -v "#" allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.vcf | cut -f 1,2 > allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_snps_1_100_K=15.txt

#A continuación, importe la información de los SNPs.

SNPb=read.table("allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_snps_1_100_K=15.txt", header=FALSE)

#Combinar los nombres de los valores atípicos con los resultados del marco de datos del bayescan.

bayescan=cbind(SNPb, bayescan)

#Cambiar el nombre de las columnas.

colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST")

#Escribe los resultados.

write.table(bayescan, "allSamples_0.90_filtred_withoutLD_bayescan_results_1_100_K=15.txt", quote=FALSE, sep="\t", row.names=FALSE)

#Cambia el valor de la columna Q_VALUE: 0 == 0.0001.

attach(bayescan)
class(bayescan$Q_VALUE)

## [1] "numeric"

bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE)
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001

#Redondea los valores.

bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4))
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4))
bayescan$ALPHA <- (round(bayescan$ALPHA, 4))
bayescan$FST <- (round(bayescan$FST, 6))

#Añada una columna para el tipo de agrupación de la selección basada en un VALOR-Q < 0,01. También puede elegir un Q-VALOR < 0,05 si quiere ser menos conservador.

bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.01,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.01,"neutral","balancing"))
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION)

#Guarde los resultados de los SNPs potencialmente bajo selección positiva (divergente) y de equilibrio (valor q < 0,01).

positive <- bayescan[bayescan$SELECTION=="diversifying",]
neutral <- bayescan[bayescan$SELECTION=="neutral",]
balancing <- bayescan[bayescan$SELECTION=="balancing",]

#Compruebe el número de SNPs que pertenecen a cada categoría.

xtabs(data=bayescan, ~SELECTION)

#SELECTION
#   balancing diversifying      neutral
#        2800          113         2610

#Escriba los resultados de los SNPs potencialmente bajo selección (valor q < 0,01).

write.table(neutral, "allSamples_0.90_1_100_K=15_neutral.txt", row.names=F, quote=F)
write.table(balancing, "allSamples_0.90_1_100_K=15_balancing.txt", row.names=F, quote=F)
write.table(positive, "allSamples_0.90_1_100_K=15_positive.txt", row.names=F, quote=F)

#Transformación Log del valor Q para crear el gráfico ggplot.

range(bayescan$Q_VALUE)

#[1] 0.0001 0.9558

bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE)

#Utilice ggplot para crear un bonito gráfico

#Crear titulos para el GRAFICO

x_title="Log(q-value)"
y_title="Fst"

#Hacer el grafico con ggplot2

library(ggplot2)

ggplot(bayescan,aes(x=LOG10_Q,y=FST)) +
  geom_point(aes(fill=SELECTION), pch=21, size=2)+
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+
  labs(x=x_title)+
  labs(y=y_title)+
  theme_classic()

#Guardar en PDF el grafico

ggsave("allSamples_0.90_bayescan.pdf", dpi=600, width=5, height=5)
dev.off()

#Utilice la función disponible en BayeScan

#También puede utilizar simplemente la función plot_R.r ya disponible en BayeScan para ver los resultados sin delimitar Q_VALUE. Primero cargue la función en R.

source ("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/R_functions/plot_R.r")

plot_bayescan("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/allSamples_0.90_withoutLD_bayescan_1_100_K=15/allSamples_0.90_filtredwithoutLD_bayescan_1_100_fst.txt")

$outliers
  [1]   15   17   48   73  174  281  288  315  320  357  378  436  461  471  486  487  579  581  603  609  673
 [22]  689  708  715  716  792  794  849  859  934  958  973  983 1008 1091 1092 1099 1107 1110 1135 1210 1219
 [43] 1230 1337 1405 1439 1460 1619 1656 1727 1766 1955 2017 2048 2094 2103 2148 2149 2160 2229 2236 2251 2304
 [64] 2315 2316 2407 2458 2513 2565 2599 2666 2679 2728 2753 2850 2868 2870 2906 2907 2908 2909 2931 2933 2934
 [85] 3009 3014 3020 3033 3099 3156 3177 3225 3282 3297 3323 3412 3512 3734 3740 3765 3794 3808 3813 3822 3823
[106] 3824 3926 3931 3996 4007 4019 4020 4126 4196 4285 4310 4403 4415 4437 4482 4502 4529 4533 4626 4646 4674
[127] 4715 4717 4718 4720 4721 4725 4729 4734 4739 4740 4750 4759 4770 4826 4836 4837 4847 4872 4898 5017 5041
[148] 5104 5105 5108 5169 5236 5237 5247 5295 5302 5338 5339 5419 5515 5516

$nb_outliers
[1] 161

Outliers_Bayescan <- positive$SNP
length (Outliers_Bayescan)
#113 outliers que se sacan de Bayescan 


#1. Descargas plink y en la misma carpeta pones tu vcf y corres la siguiente linea:

#VCF=/Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.vcf

#./plink2 --vcf $VCF --allow-extra-chr --make-bed --out /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.bed

######################################################## PCADAPT ###############################################

#Descargar la librerias

install.packages ("viridisLite")

library(pcadapt)
library(vcfR)
library(ggplot2)
library(dplyr)
library(viridis)

#Descargue los datos con la función read.pcadapt. # funciono sin el .bed se usa el vcf

setwd("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/allSamples_0.90_withoutLD_bayescan_1_100_K=15")

data <- read.pcadapt("allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.vcf", type = "vcf")

data2 <- read.vcfR("allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.vcf")

#Scanning file to determine attributes.
#File attributes:
#  meta lines: 2190
#  header_line: 2191
#  variant count: 5523
#  column count: 200
#Meta line 2190 read in.
#All meta lines processed.
#gt matrix initialized.
#Processed variant: 5523
#All variants processed
#3

data3 <- vcfR2genlight(data2)
snp <- as.data.frame(data3@loc.names)
ind <- as.data.frame(data3@ind.names)
colnames(ind) <-"INDIVIDUALS"

#Añadir información del mapa de población.

setwd("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/allSamples_0.90_withoutLD_bayescan")
poplist.names <- c(rep("AIRE",1), rep("ESCURT", 1), rep("ADRG", 1), rep("MOLTONA", 1), rep("ESCURT",22), rep("MOLTONA",10), rep("GUARDIA",15), rep("AIRE",10), rep("COLOM",10), rep("PRCV",10), rep("ADRG",9), rep("ROVELLS",10), rep("REI",1), rep("AIRE", 1), rep("REI", 9), rep("PORROS", 1), rep("AIRE",1), rep("PORROS",10), rep("AIRE",1), rep("FORADADA",10), rep("AIRE",1), rep("ESCLTS",10), rep("AIRE",1), rep("DRAGONERA",4), rep("CABRERAH",6), rep("AIRE",1), rep("CABRERAH",6), rep("COLOM",4), rep("AIRE",1), rep("COLOM",6), rep("CABREBRAI",4), rep("AIRE",1), rep("CABRERAI",4), rep("COLOMER",4), rep("AIRE",4))
print(poplist.names)


#Primero ejecute pcadapt con un gran número K de componentes principales, por ejemplo K = 20

data_pcadapt_trial <- pcadapt(data, K = 20, min.maf = 0.01, ploidy=2)

#Haga un screeplot y un score plot para determinar el número ideal de PC. El 'scree plot' muestra en orden decreciente el porcentaje de varianza explicado por cada PC.

plot(data_pcadapt_trial, option = "screeplot", col="blue", snp.info = NULL,
     plt.pkg = "ggplot")

#Compruebe el gráfico para seleccionar el número óptimo K de componentes principales. En primer lugar, compruebe PC1 y PC2.

plot(data_pcadapt_trial, option = "scores", pop=poplist.names,
     snp.info = NULL, plt.pkg = "ggplot")

# Vamos a chequear el PC3 Y PC4

plot(data_pcadapt_trial, option = "scores", pop=poplist.names,
     i = 3, j = 4, snp.info = NULL, plt.pkg = "ggplot")

#Vamos a chequear PC5 Y PC6

plot(data_pcadapt_trial, option = "scores", pop=poplist.names,
     i = 5, j = 6, snp.info = NULL, plt.pkg = "ggplot")

#Calcular la prueba estadística

#El estadístico de prueba para detectar SNPs atípicos es la distancia de Mahalanobis, que es un enfoque multidimensional que mide lo distante que está un punto de la media.

#Para nuestro conjunto de datos, ejecutaremos pcadapt con K = 15 (ya que Se recomienda mantener los PC que corresponden a los valores propios a la izquierda de la línea recta (regla de Cattell).

#Por defecto, el parámetro min.maf está fijado en el 5%; aquí lo cambiamos al 1%. Los valores p de los SNPs con una frecuencia alélica menor que el umbral no se calculan (se devuelve NA).

data_pcadapt <- pcadapt(data, K = 15, min.maf = 0.01)

#Obtener una lista de valores P.

snps_pvalues <- cbind(snp, data_pcadapt$pvalues)
snps_pvalues_no_na <- na.omit(snps_pvalues)
write.table(snps_pvalues, "Allsamples_0.90_withoutLD_Pvalues_1_100_K=15_.txt", sep="\t", quote=FALSE)
write.table(snps_pvalues_no_na, "Allsamples_0.90_withoutLD_Pvalues_SINNA_1_100_K=15_.txt", sep="\t", quote=FALSE)

#Herramientas gráficas: el usuario también puede comprobar la distribución uniforme esperada de los valores p mediante un gráfico Q-Q

plot(data_pcadapt, option = "qqplot", col, snp.info = NULL, plt.pkg = "ggplot")

#Aquí el gráfico Q-Q confirma que la mayoría de los valores p siguen la distribución uniforme esperada.

#Compruebe la distribución de estos valores p.

hist(data_pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

#Un histograma de los valores p confirma que la mayoría de los valores p siguen una distribución uniforme. El exceso de valores P pequeños indica la presencia de valores atípicos.

#Visualizar la distribución de los valores p

quantile(snps_pvalues_no_na$`data_pcadapt$pvalues`, probs = c(0.01, 0.99))

#1%          99%
#2.323840e-15 9.999674e-01

#Una opcion para hacer Cutt off (Pero no se deja para el paper, se usara qvalue < 0.01)
#Obtenga sólo los marcadores que muestran valores p extremos: el 1% superior. #OJOOOOOOOO CAMBIAR EL NUMERO
#OJOOOOOOOOO

#top_1percent <- subset(snps_pvalues_no_na, snps_pvalues_no_na$`data_pcadapt$pvalues` <= 2.323840e-15)
#colnames(top_1percent) <- c("LOCUS","PVALUE")
#write.table(top_1percent, "allSamples_0.90_withoutLD_outliers_pcadapt_1_100_K=15.txt", sep="\t", quote=FALSE, row.names = FALSE)
#length (top_1percent$LOCUS)
#50

##################### q-values to Cutoff  #######################

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("qvalue")

library("qvalue")
snps_pvalues_no_na$qval <- qvalue(snps_pvalues_no_na$`data_pcadapt$pvalues`)$qvalues
alpha <- 0.01
outliers_qvalue<- snps_pvalues_no_na$`data3@loc.names`[which(snps_pvalues_no_na$qval < alpha)]
length(outliers_qvalue)
#442
write.table(outliers_qvalue, "allSamples_0.90_withoutLD_outliers_pcadapt_1_100_K=15.txt", sep="\t", quote=FALSE, row.names = FALSE)

################## Otras opciones de Cutt off ##################################
###################Benjamini-Hochberg Procedure###############################

#padj <- p.adjust(snps_pvalues_no_na$pvalues,method="BH")
#alpha <- 0.05
#outliers_benjamini <- which(padj < alpha)
#length(outliers_benjamini)
#650

########### Bonferroni correction ############

#padj_1 <- p.adjust(data_pcadapt$pvalues,method="bonferroni")
#alpha <- 0.1
#outliers_bonferroni <- which(padj_1 < alpha)
#length(outliers_bonferroni)
#242


####################################################### VENN DIAGRAM ###################################

#Cargar las librerias

#install.packages("VennDiagram")
library(VennDiagram)

#Guardar los nombres de SNP atípicos de Bayescan en un vector.
bayescan_outliers <- positive$SNP

#Compruebe el número de valores atípicos encontrados por BayeScan.

length(bayescan_outliers)
#[1] 113

#Guarda los nombres de pcadapt en un vector.

pcadapt_outliers <- outliers_qvalue

#Comprueba el número de valores atípicos encontrados por pcadapt.

length(pcadapt_outliers)
#[1] 442

#Crea un diagrama de VENN

venn.diagram(
x = list(bayescan_outliers,pcadapt_outliers),
category.names = c("Bayescan" , "Pcadapt"),
filename = 'Venn_diagramm_outliers.png',
output=TRUE,
imagetype="png" ,
height = 400 ,
width = 400,
resolution = 300,
compression = "lzw",
cat.cex = 0.6,
cat.pos = c(-5, 5))

#El grafico se guarda en la carpeta directamente

#Diagrama de Venn.
#Encuentre los valores atípicos comunes a ambos BayeScan y pcadapt.

outliers_comunes <- intersect(bayescan_outliers,pcadapt_outliers)
write.table(outliers_comunes, "outliers_comunes_allSamples_0.90_1_100_K=15.txt", row.names=F, quote=F)

[1] "rPodLil1.1_Chr1_2934019"   "rPodLil1.1_Chr1_46895223" 
[3] "rPodLil1.1_Chr1_95267778"  "rPodLil1.1_Chr1_117517707"
[5] "rPodLil1.1_Chr1_119179273" "rPodLil1.1_Chr2_24907248" 
[7] "rPodLil1.1_Chr2_27008399"  "rPodLil1.1_Chr2_28593619" 
[9] "rPodLil1.1_Chr2_42088672"  "rPodLil1.1_Chr2_42095785" 
[11] "rPodLil1.1_Chr2_59052307"  "rPodLil1.1_Chr2_84098004" 
[13] "rPodLil1.1_Chr2_113395694" "rPodLil1.1_Chr2_113395780"
[15] "rPodLil1.1_Chr2_113677261" "rPodLil1.1_Chr2_115455674"
[17] "rPodLil1.1_Chr3_3617630"   "rPodLil1.1_Chr4_66903019" 
[19] "rPodLil1.1_Chr4_78265967"  "rPodLil1.1_Chr5_46265915" 
[21] "rPodLil1.1_Chr5_84130007"  "rPodLil1.1_Chr6_1795447"  
[23] "rPodLil1.1_Chr6_59810146"  "rPodLil1.1_Chr6_90332840" 
[25] "rPodLil1.1_Chr7_35098419"  "rPodLil1.1_Chr7_42515732" 
[27] "rPodLil1.1_Chr7_88655001"  "rPodLil1.1_Chr7_88655031" 
[29] "rPodLil1.1_Chr7_88655044"  "rPodLil1.1_Chr7_88956143" 
[31] "rPodLil1.1_Chr8_5488936"   "rPodLil1.1_Chr8_26520893" 
[33] "rPodLil1.1_Chr8_27743549"  "rPodLil1.1_Chr8_28398888" 
[35] "rPodLil1.1_Chr8_69804324"  "rPodLil1.1_Chr10_55226525"
[37] "rPodLil1.1_Chr11_1468675"  "rPodLil1.1_Chr11_39073477"
[39] "rPodLil1.1_Chr11_59491095" "rPodLil1.1_Chr12_862930"  
[41] "rPodLil1.1_Chr12_937545"   "rPodLil1.1_Chr13_4577558" 
[43] "rPodLil1.1_Chr13_12207764" "rPodLil1.1_Chr14_50270897"
[45] "rPodLil1.1_Chr14_51265174" "rPodLil1.1_Chr14_51265175"
[47] "rPodLil1.1_Chr14_51265204" "rPodLil1.1_Chr14_52332112"
[49] "rPodLil1.1_ChrZ_3302964"   "rPodLil1.1_ChrZ_3640603"  
[51] "rPodLil1.1_ChrZ_3727349"   "rPodLil1.1_ChrZ_8962162"  
[53] "rPodLil1.1_ChrZ_10885052"  "rPodLil1.1_ChrZ_13397106" 
[55] "rPodLil1.1_ChrZ_35140541"  "rPodLil1.1_ChrZ_35140565" 
[57] "rPodLil1.1_Chr15_43227153" "rPodLil1.1_Chr15_43987288"
[59] "rPodLil1.1_Chr16_19655125" "rPodLil1.1_Chr16_31685302"

length(intersect(bayescan_outliers,pcadapt_outliers))
#60

#Como se puede ver, el subconjunto del número compartido de candidatos es bajo, lo que sugiere que los métodos de exploración del genoma son muy variables.
#Este hallazgo está de acuerdo con Dalongeville et al. 2018, quienes mostraron que la identificación y el número de los valores atípicos son drásticamente diferentes dependiendo del método utilizado.

#Encuentre los valores atípicos únicos a BayeScan.

setdiff(bayescan_outliers,pcadapt_outliers)

[1] "rPodLil1.1_Chr1_3198243"   "rPodLil1.1_Chr1_75194630" 
[3] "rPodLil1.1_Chr1_78815279"  "rPodLil1.1_Chr1_109795225"
[5] "rPodLil1.1_Chr1_122505450" "rPodLil1.1_Chr1_122505496"
[7] "rPodLil1.1_Chr2_373437"    "rPodLil1.1_Chr2_1480770"  
[9] "rPodLil1.1_Chr2_9626526"   "rPodLil1.1_Chr2_78426729" 
[11] "rPodLil1.1_Chr2_94422121"  "rPodLil1.1_Chr3_30325923" 
[13] "rPodLil1.1_Chr3_33351491"  "rPodLil1.1_Chr3_36947906" 
[15] "rPodLil1.1_Chr3_105358002" "rPodLil1.1_Chr4_31108867" 
[17] "rPodLil1.1_Chr5_28312696"  "rPodLil1.1_Chr5_54826928" 
[19] "rPodLil1.1_Chr5_69862185"  "rPodLil1.1_Chr5_82969696" 
[21] "rPodLil1.1_Chr6_25364734"  "rPodLil1.1_Chr6_25364773" 
[23] "rPodLil1.1_Chr6_47059684"  "rPodLil1.1_Chr6_70321862" 
[25] "rPodLil1.1_Chr7_20250812"  "rPodLil1.1_Chr7_81173402" 
[27] "rPodLil1.1_Chr7_82758419"  "rPodLil1.1_Chr8_5488947"  
[29] "rPodLil1.1_Chr8_5488953"   "rPodLil1.1_Chr8_59036780" 
[31] "rPodLil1.1_Chr8_78143746"  "rPodLil1.1_Chr8_81521452" 
[33] "rPodLil1.1_Chr9_2239317"   "rPodLil1.1_Chr10_69972391"
[35] "rPodLil1.1_Chr11_1468708"  "rPodLil1.1_Chr11_1468711" 
[37] "rPodLil1.1_Chr11_58136094" "rPodLil1.1_Chr12_29780262"
[39] "rPodLil1.1_Chr12_47862849" "rPodLil1.1_Chr13_42731388"
[41] "rPodLil1.1_Chr14_7801041"  "rPodLil1.1_Chr14_16581151"
[43] "rPodLil1.1_Chr14_18249460" "rPodLil1.1_Chr14_43543747"
[45] "rPodLil1.1_Chr14_51265201" "rPodLil1.1_ChrZ_39439943" 
[47] "rPodLil1.1_ChrZ_47120167"  "rPodLil1.1_Chr15_21475702"
[49] "rPodLil1.1_Chr15_26306341" "rPodLil1.1_Chr17_1846576" 
[51] "rPodLil1.1_Chr16_2388980"  "rPodLil1.1_Chr16_16086701"
[53] "rPodLil1.1_Chr16_19655100"

#Ver el numero de Outliers enconrtrado por Bayescan

length(setdiff(bayescan_outliers,pcadapt_outliers))
#[1] 53

#Encuentre los valores atípicos únicos para pcadaptarse.

setdiff(pcadapt_outliers,bayescan_outliers)

[1] "rPodLil1.1_Chr1_3953711"   "rPodLil1.1_Chr1_3953749"  
[3] "rPodLil1.1_Chr1_3953753"   "rPodLil1.1_Chr1_4284904"  
[5] "rPodLil1.1_Chr1_7275620"   "rPodLil1.1_Chr1_31366687" 
[7] "rPodLil1.1_Chr1_44135202"  "rPodLil1.1_Chr1_44470512" 
[9] "rPodLil1.1_Chr1_44470531"  "rPodLil1.1_Chr1_52994644" 
[11] "rPodLil1.1_Chr1_58076355"  "rPodLil1.1_Chr1_58558759" 
[13] "rPodLil1.1_Chr1_61403110"  "rPodLil1.1_Chr1_71720582" 
[15] "rPodLil1.1_Chr1_75905596"  "rPodLil1.1_Chr1_78815237" 
[17] "rPodLil1.1_Chr1_87965406"  "rPodLil1.1_Chr1_89874910" 
[19] "rPodLil1.1_Chr1_95368212"  "rPodLil1.1_Chr1_97258310" 
[21] "rPodLil1.1_Chr1_99103946"  "rPodLil1.1_Chr1_107469876"
[23] "rPodLil1.1_Chr1_119091652" "rPodLil1.1_Chr1_123429581"
[25] "rPodLil1.1_Chr1_125813798" "rPodLil1.1_Chr1_126232054"
[27] "rPodLil1.1_Chr1_127309209" "rPodLil1.1_Chr1_134802137"
[29] "rPodLil1.1_Chr1_134947894" "rPodLil1.1_Chr1_136748232"
[31] "rPodLil1.1_Chr2_22550367"  "rPodLil1.1_Chr2_25087310" 
[33] "rPodLil1.1_Chr2_27438773"  "rPodLil1.1_Chr2_29225230" 
[35] "rPodLil1.1_Chr2_29396781"  "rPodLil1.1_Chr2_34417212" 
[37] "rPodLil1.1_Chr2_38226719"  "rPodLil1.1_Chr2_40437107" 
[39] "rPodLil1.1_Chr2_40927743"  "rPodLil1.1_Chr2_40927795" 
[41] "rPodLil1.1_Chr2_41493724"  "rPodLil1.1_Chr2_42088704" 
[43] "rPodLil1.1_Chr2_48189787"  "rPodLil1.1_Chr2_48189801" 
[45] "rPodLil1.1_Chr2_48189845"  "rPodLil1.1_Chr2_51622620" 
[47] "rPodLil1.1_Chr2_54239377"  "rPodLil1.1_Chr2_55970079" 
[49] "rPodLil1.1_Chr2_56946110"  "rPodLil1.1_Chr2_57819665" 
[51] "rPodLil1.1_Chr2_60641496"  "rPodLil1.1_Chr2_70684018" 
[53] "rPodLil1.1_Chr2_72827927"  "rPodLil1.1_Chr2_84790682" 
[55] "rPodLil1.1_Chr2_85831218"  "rPodLil1.1_Chr2_87001885" 
[57] "rPodLil1.1_Chr2_89916085"  "rPodLil1.1_Chr2_98047062" 
[59] "rPodLil1.1_Chr2_102572719" "rPodLil1.1_Chr2_106223568"
[61] "rPodLil1.1_Chr2_109912296" "rPodLil1.1_Chr2_110736721"
[63] "rPodLil1.1_Chr2_111228264" "rPodLil1.1_Chr2_111813741"
[65] "rPodLil1.1_Chr2_113614150" "rPodLil1.1_Chr2_114081815"
[67] "rPodLil1.1_Chr2_114081842" "rPodLil1.1_Chr2_115308818"
[69] "rPodLil1.1_Chr2_117754647" "rPodLil1.1_Chr2_120859984"
[71] "rPodLil1.1_Chr2_120880260" "rPodLil1.1_Chr2_122344747"
[73] "rPodLil1.1_Chr2_122344749" "rPodLil1.1_Chr2_122991008"
[75] "rPodLil1.1_Chr2_122991062" "rPodLil1.1_Chr3_445294"   
[77] "rPodLil1.1_Chr3_938013"    "rPodLil1.1_Chr3_22818955" 
[79] "rPodLil1.1_Chr3_29820163"  "rPodLil1.1_Chr3_38164532" 
[81] "rPodLil1.1_Chr3_49700906"  "rPodLil1.1_Chr3_54547934" 
[83] "rPodLil1.1_Chr3_68533812"  "rPodLil1.1_Chr3_68788039" 
[85] "rPodLil1.1_Chr3_69930547"  "rPodLil1.1_Chr3_70136544" 
[87] "rPodLil1.1_Chr3_73473673"  "rPodLil1.1_Chr3_77388789" 
[89] "rPodLil1.1_Chr3_81012436"  "rPodLil1.1_Chr3_81012455" 
[91] "rPodLil1.1_Chr3_81012479"  "rPodLil1.1_Chr3_95003289" 
[93] "rPodLil1.1_Chr3_95003319"  "rPodLil1.1_Chr3_97717445" 
[95] "rPodLil1.1_Chr3_98028234"  "rPodLil1.1_Chr3_100360121"
[97] "rPodLil1.1_Chr3_102141064" "rPodLil1.1_Chr3_105283519"
[99] "rPodLil1.1_Chr3_105358004" "rPodLil1.1_Chr3_108640658"
[101] "rPodLil1.1_Chr3_109541494" "rPodLil1.1_Chr3_113472345"
[103] "rPodLil1.1_Chr3_113666086" "rPodLil1.1_Chr3_113678600"
[105] "rPodLil1.1_Chr3_118246051" "rPodLil1.1_Chr3_119155518"
[107] "rPodLil1.1_Chr4_2799338"   "rPodLil1.1_Chr4_33700304" 
[109] "rPodLil1.1_Chr4_37291859"  "rPodLil1.1_Chr4_42004939" 
[111] "rPodLil1.1_Chr4_49372586"  "rPodLil1.1_Chr4_57635480" 
[113] "rPodLil1.1_Chr4_61614452"  "rPodLil1.1_Chr4_75353229" 
[115] "rPodLil1.1_Chr4_75353255"  "rPodLil1.1_Chr4_75407352" 
[117] "rPodLil1.1_Chr4_77872702"  "rPodLil1.1_Chr4_81738033" 
[119] "rPodLil1.1_Chr4_85915050"  "rPodLil1.1_Chr4_85915104" 
[121] "rPodLil1.1_Chr4_86706839"  "rPodLil1.1_Chr4_93898698" 
[123] "rPodLil1.1_Chr4_102526255" "rPodLil1.1_Chr5_4305688"  
[125] "rPodLil1.1_Chr5_12403679"  "rPodLil1.1_Chr5_18905612" 
[127] "rPodLil1.1_Chr5_19497405"  "rPodLil1.1_Chr5_22866818" 
[129] "rPodLil1.1_Chr5_24603327"  "rPodLil1.1_Chr5_25182047" 
[131] "rPodLil1.1_Chr5_27633596"  "rPodLil1.1_Chr5_29498015" 
[133] "rPodLil1.1_Chr5_29654553"  "rPodLil1.1_Chr5_35927836" 
[135] "rPodLil1.1_Chr5_55676888"  "rPodLil1.1_Chr5_56705486" 
[137] "rPodLil1.1_Chr5_68981555"  "rPodLil1.1_Chr5_69862168" 
[139] "rPodLil1.1_Chr5_78409513"  "rPodLil1.1_Chr5_80477184" 
[141] "rPodLil1.1_Chr5_82351577"  "rPodLil1.1_Chr5_91035926" 
[143] "rPodLil1.1_Chr5_91182507"  "rPodLil1.1_Chr5_91182536" 
[145] "rPodLil1.1_Chr5_96058849"  "rPodLil1.1_Chr6_1870499"  
[147] "rPodLil1.1_Chr6_1870510"   "rPodLil1.1_Chr6_14751355" 
[149] "rPodLil1.1_Chr6_14967056"  "rPodLil1.1_Chr6_14967059" 
[151] "rPodLil1.1_Chr6_19435223"  "rPodLil1.1_Chr6_24482651" 
[153] "rPodLil1.1_Chr6_29713507"  "rPodLil1.1_Chr6_44209253" 
[155] "rPodLil1.1_Chr6_50599339"  "rPodLil1.1_Chr6_50599353" 
[157] "rPodLil1.1_Chr6_53067906"  "rPodLil1.1_Chr6_54340681" 
[159] "rPodLil1.1_Chr6_55421868"  "rPodLil1.1_Chr6_59360390" 
[161] "rPodLil1.1_Chr6_63589717"  "rPodLil1.1_Chr6_64639249" 
[163] "rPodLil1.1_Chr6_75271970"  "rPodLil1.1_Chr6_78696740" 
[165] "rPodLil1.1_Chr6_80203961"  "rPodLil1.1_Chr7_3902043"  
[167] "rPodLil1.1_Chr7_5480271"   "rPodLil1.1_Chr7_5798565"  
[169] "rPodLil1.1_Chr7_9285978"   "rPodLil1.1_Chr7_10047863" 
[171] "rPodLil1.1_Chr7_15228236"  "rPodLil1.1_Chr7_22738031" 
[173] "rPodLil1.1_Chr7_25760385"  "rPodLil1.1_Chr7_26973889" 
[175] "rPodLil1.1_Chr7_29229822"  "rPodLil1.1_Chr7_29682087" 
[177] "rPodLil1.1_Chr7_45158707"  "rPodLil1.1_Chr7_57768346" 
[179] "rPodLil1.1_Chr7_63254920"  "rPodLil1.1_Chr7_74939318" 
[181] "rPodLil1.1_Chr7_78882285"  "rPodLil1.1_Chr7_80021068" 
[183] "rPodLil1.1_Chr7_81173352"  "rPodLil1.1_Chr7_82263253" 
[185] "rPodLil1.1_Chr7_85455243"  "rPodLil1.1_Chr7_86717694" 
[187] "rPodLil1.1_Chr7_87930217"  "rPodLil1.1_Chr8_316833"   
[189] "rPodLil1.1_Chr8_5488987"   "rPodLil1.1_Chr8_11384983" 
[191] "rPodLil1.1_Chr8_11491300"  "rPodLil1.1_Chr8_11956960" 
[193] "rPodLil1.1_Chr8_12258640"  "rPodLil1.1_Chr8_25026159" 
[195] "rPodLil1.1_Chr8_25508299"  "rPodLil1.1_Chr8_28084370" 
[197] "rPodLil1.1_Chr8_35329907"  "rPodLil1.1_Chr8_36427689" 
[199] "rPodLil1.1_Chr8_36569209"  "rPodLil1.1_Chr8_44773616" 
[201] "rPodLil1.1_Chr8_45584398"  "rPodLil1.1_Chr8_45966861" 
[203] "rPodLil1.1_Chr8_47477215"  "rPodLil1.1_Chr8_47477446" 
[205] "rPodLil1.1_Chr8_48629135"  "rPodLil1.1_Chr8_60362332" 
[207] "rPodLil1.1_Chr8_65153249"  "rPodLil1.1_Chr8_66139439" 
[209] "rPodLil1.1_Chr8_70164825"  "rPodLil1.1_Chr8_71843103" 
[211] "rPodLil1.1_Chr8_73345930"  "rPodLil1.1_Chr8_73433128" 
[213] "rPodLil1.1_Chr8_73814065"  "rPodLil1.1_Chr8_73923326" 
[215] "rPodLil1.1_Chr8_76870977"  "rPodLil1.1_Chr8_81587371" 
[217] "rPodLil1.1_Chr8_82841842"  "rPodLil1.1_Chr9_2376332"  
[219] "rPodLil1.1_Chr9_5113219"   "rPodLil1.1_Chr9_8115813"  
[221] "rPodLil1.1_Chr9_14503863"  "rPodLil1.1_Chr9_19560609" 
[223] "rPodLil1.1_Chr9_24004247"  "rPodLil1.1_Chr9_24027092" 
[225] "rPodLil1.1_Chr9_24056004"  "rPodLil1.1_Chr9_24346306" 
[227] "rPodLil1.1_Chr9_29145237"  "rPodLil1.1_Chr9_29145263" 
[229] "rPodLil1.1_Chr9_41751849"  "rPodLil1.1_Chr9_45150546" 
[231] "rPodLil1.1_Chr9_49740089"  "rPodLil1.1_Chr9_54101520" 
[233] "rPodLil1.1_Chr9_55376629"  "rPodLil1.1_Chr9_56641259" 
[235] "rPodLil1.1_Chr9_57134091"  "rPodLil1.1_Chr9_57134125" 
[237] "rPodLil1.1_Chr9_58900168"  "rPodLil1.1_Chr9_60381040" 
[239] "rPodLil1.1_Chr9_62487656"  "rPodLil1.1_Chr9_69236611" 
[241] "rPodLil1.1_Chr9_69236686"  "rPodLil1.1_Chr9_74308700" 
[243] "rPodLil1.1_Chr9_74308727"  "rPodLil1.1_Chr9_74308751" 
[245] "rPodLil1.1_Chr10_444099"   "rPodLil1.1_Chr10_6672186" 
[247] "rPodLil1.1_Chr10_53696153" "rPodLil1.1_Chr10_56345717"
[249] "rPodLil1.1_Chr10_70976079" "rPodLil1.1_Chr11_27706553"
[251] "rPodLil1.1_Chr11_36935087" "rPodLil1.1_Chr11_37359286"
[253] "rPodLil1.1_Chr11_39780652" "rPodLil1.1_Chr11_39780664"
[255] "rPodLil1.1_Chr11_56743314" "rPodLil1.1_Chr12_4185400" 
[257] "rPodLil1.1_Chr12_14784108" "rPodLil1.1_Chr12_17721543"
[259] "rPodLil1.1_Chr12_23537490" "rPodLil1.1_Chr12_29780264"
[261] "rPodLil1.1_Chr12_33001656" "rPodLil1.1_Chr12_37918224"
[263] "rPodLil1.1_Chr12_40512688" "rPodLil1.1_Chr12_41129150"
[265] "rPodLil1.1_Chr12_41816290" "rPodLil1.1_Chr12_45178767"
[267] "rPodLil1.1_Chr12_48830361" "rPodLil1.1_Chr12_50270936"
[269] "rPodLil1.1_Chr12_50829051" "rPodLil1.1_Chr12_50829064"
[271] "rPodLil1.1_Chr13_496470"   "rPodLil1.1_Chr13_1200400" 
[273] "rPodLil1.1_Chr13_9477984"  "rPodLil1.1_Chr13_13230157"
[275] "rPodLil1.1_Chr13_13375862" "rPodLil1.1_Chr13_13396663"
[277] "rPodLil1.1_Chr13_13396689" "rPodLil1.1_Chr13_16874356"
[279] "rPodLil1.1_Chr13_16970051" "rPodLil1.1_Chr13_17254725"
[281] "rPodLil1.1_Chr13_18320462" "rPodLil1.1_Chr13_22827420"
[283] "rPodLil1.1_Chr13_24596268" "rPodLil1.1_Chr13_45203598"
[285] "rPodLil1.1_Chr14_1705777"  "rPodLil1.1_Chr14_1773054" 
[287] "rPodLil1.1_Chr14_2679404"  "rPodLil1.1_Chr14_5308954" 
[289] "rPodLil1.1_Chr14_9285326"  "rPodLil1.1_Chr14_9449883" 
[291] "rPodLil1.1_Chr14_9449900"  "rPodLil1.1_Chr14_9562564" 
[293] "rPodLil1.1_Chr14_9562623"  "rPodLil1.1_Chr14_10022818"
[295] "rPodLil1.1_Chr14_13863860" "rPodLil1.1_Chr14_22506358"
[297] "rPodLil1.1_Chr14_23882523" "rPodLil1.1_Chr14_24511253"
[299] "rPodLil1.1_Chr14_24533055" "rPodLil1.1_Chr14_24786146"
[301] "rPodLil1.1_Chr14_24786184" "rPodLil1.1_Chr14_25141883"
[303] "rPodLil1.1_Chr14_25141911" "rPodLil1.1_Chr14_25213773"
[305] "rPodLil1.1_Chr14_35657702" "rPodLil1.1_Chr14_36306094"
[307] "rPodLil1.1_Chr14_41325965" "rPodLil1.1_Chr14_45435100"
[309] "rPodLil1.1_Chr14_46734306" "rPodLil1.1_Chr14_49172258"
[311] "rPodLil1.1_Chr14_49172270" "rPodLil1.1_Chr14_51265197"
[313] "rPodLil1.1_Chr14_51994559" "rPodLil1.1_Chr14_52331964"
[315] "rPodLil1.1_ChrZ_6834001"   "rPodLil1.1_ChrZ_10885020" 
[317] "rPodLil1.1_ChrZ_12322318"  "rPodLil1.1_ChrZ_12960668" 
[319] "rPodLil1.1_ChrZ_15108107"  "rPodLil1.1_ChrZ_21784230" 
[321] "rPodLil1.1_ChrZ_26929334"  "rPodLil1.1_ChrZ_32486707" 
[323] "rPodLil1.1_ChrZ_32486744"  "rPodLil1.1_ChrZ_35530766" 
[325] "rPodLil1.1_ChrZ_41763205"  "rPodLil1.1_ChrZ_42458098" 
[327] "rPodLil1.1_ChrZ_42458204"  "rPodLil1.1_ChrZ_44573524" 
[329] "rPodLil1.1_ChrZ_47506320"  "rPodLil1.1_ChrZ_48555853" 
[331] "rPodLil1.1_Chr15_568690"   "rPodLil1.1_Chr15_2262395" 
[333] "rPodLil1.1_Chr15_5623619"  "rPodLil1.1_Chr15_8688673" 
[335] "rPodLil1.1_Chr15_8688686"  "rPodLil1.1_Chr15_8737450" 
[337] "rPodLil1.1_Chr15_12448676" "rPodLil1.1_Chr15_19638436"
[339] "rPodLil1.1_Chr15_21215115" "rPodLil1.1_Chr15_27637247"
[341] "rPodLil1.1_Chr15_27837231" "rPodLil1.1_Chr15_29110728"
[343] "rPodLil1.1_Chr15_32469686" "rPodLil1.1_Chr15_32469711"
[345] "rPodLil1.1_Chr15_34123591" "rPodLil1.1_Chr15_38964861"
[347] "rPodLil1.1_Chr15_41465179" "rPodLil1.1_Chr17_3760798" 
[349] "rPodLil1.1_Chr17_4181736"  "rPodLil1.1_Chr17_21659547"
[351] "rPodLil1.1_Chr17_37126830" "rPodLil1.1_Chr17_39858061"
[353] "rPodLil1.1_Chr16_1597075"  "rPodLil1.1_Chr16_1727130" 
[355] "rPodLil1.1_Chr16_1727175"  "rPodLil1.1_Chr16_2684644" 
[357] "rPodLil1.1_Chr16_3593348"  "rPodLil1.1_Chr16_6944562" 
[359] "rPodLil1.1_Chr16_11864343" "rPodLil1.1_Chr16_13483744"
[361] "rPodLil1.1_Chr16_16738550" "rPodLil1.1_Chr16_16828993"
[363] "rPodLil1.1_Chr16_19036845" "rPodLil1.1_Chr16_20268420"
[365] "rPodLil1.1_Chr16_20322592" "rPodLil1.1_Chr16_22342749"
[367] "rPodLil1.1_Chr16_23107320" "rPodLil1.1_Chr16_26784242"
[369] "rPodLil1.1_Chr16_29210634" "rPodLil1.1_Chr16_30271301"
[371] "rPodLil1.1_Chr16_30396100" "rPodLil1.1_Chr16_30578931"
[373] "rPodLil1.1_Chr16_30675603" "rPodLil1.1_Chr16_30675647"
[375] "rPodLil1.1_Chr16_30885887" "rPodLil1.1_Chr16_40557972"
[377] "rPodLil1.1_Chr16_40557973" "rPodLil1.1_Chr18_875270"  
[379] "rPodLil1.1_Chr18_3061384"  "rPodLil1.1_Chr18_6930344" 
[381] "rPodLil1.1_Chr18_7522829"  "rPodLil1.1_Chr18_7522873" 

#Ver el numero de Outliers enconrtrado por Bayescan
length(setdiff(pcadapt_outliers,bayescan_outliers))
#[1] 382

#IR A VCFTOOLS Y EXTRAER LOS OUTLIERS DE BayeScan
# Hacer un archivo con los outrliers de bayescan con la informacion del cromosoma un TAB y luego la posicion :  outliers_comunes_allSamples_0.90_1_100_K=15_2.txt
#En vcftools la opcion positions me permite extraer los snps que me interesan

cd /Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/allSamples_0.90_withoutLD_bayescan_1_100_K=15

vcftools --vcf allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.vcf --positions outliers_comunes_allSamples_0.90_1_100_K=15_2.txt --recode --recode-INFO-all --out outliers_comunes_bayescan_pcadapt_allSamples_0.90_1_100_K=15

grep -v # outliers_comunes_bayescan_pcadapt_allSamples_0.90_1_100_K=15.recode.vcf | wc -l
