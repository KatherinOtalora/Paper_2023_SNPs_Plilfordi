############### R PIPELINE SNPS ############

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
install.packages("hierfstat") #No se pudo instalar en R studio
install.packages("devtools")
install_github("green-striped-gecko/dartR") #No se pudo instalar en Rstudio

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
install.packages("SNPRelate")
install.packages("ggplot2")
install.packages("factoextra")

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")
install.packages("tidyr")

library(devtools)
install_github("mijangos81/ggtern", force = TRUE)
gl.install.vanilla.dartR()
install.packages("dartR")

install.packages("terra")
library("dartR")

install.packages("MASS") 
install.packages("reshape2") 
install.packages("reshape") 

#Lllamar las librerias necesarias para el analisis de SNPs

library(vcfR)
library(poppr)
library(ape)
library(pheatmap)
library(RColorBrewer)
library(pillar)
library(adegenet)
library(dartR)
library("dartR")
library(adegenet)
library(hierfstat)
library(devtools)
library(pegas)
library(ggplot2)
library(factoextra)
library(tidyr)
library(SNPRelate)
library(MASS) 
library(reshape2) 
library(reshape) 

#Nos ubicamos en la carpeta local desde donde vamos a trabajar, hay que asegúrese de que está en la carpeta correcta con los archivos descargados disponibles.

setwd("/Users/katherin_otalora/Documents/R_analysis/GATK_2022/allSamples/allSamples_0.90")

#El output de LD from PLINK tiene los individuos con doble nombre y hay que dejarlos con su nombre individual, se puede cambiar con texwangler esto es por el --double-id - told plink to duplicate the id of our samples

# Asignamos  nuestros datos .vcf filtrados por LD a una variable nueva con nombre podarcislilfordi.VCF

podarcislilfordi.VCF <- read.vcfR("allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_1.vcf")

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
gl.podarcislilfordi
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

########################## Iniciamos calculando el FST Pairwaise with Bottstrap p-value para nuestra muestra de P.lilfordi ############################

#Utilizaremos la funcion gl.fst.pop en donde hay que poner nuestro objeto genlight, nboots representan número de bootstraps a realizar entre loci para generar intervalos de confianza y valores p
#Tambien tenemos el percent que significa el percentil para calcular el intervalo de confianza alrededor de, por defalut = 95.
# La funcion nclusters representan el número de núcleos del procesador a utilizar durante los cálculos.

Pairwaise_Fst_Plilfordi <-gl.fst.pop(gl.podarcislilfordi, nboots = 1000, percent = 95, nclusters = 1)
Pairwaise_Fst_Plilfordi

########################### HEATMAP Fst ###############################
# El HetMap es una técnica de visualización de datos que muestra la magnitud de un fenómeno como color en dos dimensiones. La variación en el color puede ser por matiz o intensidad, dando señales visuales obvias al lector sobre cómo el fenómeno se agrupa o varía en el espacio.
# En este caso la vareiable que vamos a utilizar para visualizarla es el FST o diferenciacion genetica entre poblaciones de P.lilfordi

#Previo al HeatMap hay que Hacer la matrix en excel poniendo la misma forma de la tabla de FST pero totalmente rellena (Matriz espejo), guardarla desde excel en formarto .txt y luego Importanla en R
#Al guardar el .txt recordar quitarle la primera linea que dice $Fsts

setwd("/Users/katherin_otalora/Documents/R_analysis/GATK_2022/allSamples/allSamples_0.90")
df <- read.table ("FST_B10000_allSamples_GATK_22_05_23.txt", header = TRUE)
#Encurt y Porros cerca 
df <- read.table ("1_FST_B10000_allSamples_GATK_22_05_23.txt", header = TRUE)
df
df_matrix <-as.matrix(df, rownames=TRUE, rownames.value=NULL)
df_matrix
# Con dendrograma
pheatmap(df_matrix, color = hcl.colors(50, "BluYl"), display_numbers = TRUE, number_color = "black", fontsize_number = 8, cluster_rows = FALSE,)
#Sin dendograma
pheatmap(df_matrix, color = hcl.colors(50, "BluYl"), display_numbers = TRUE, number_color = "black", fontsize_number = 8, cluster_rows = FALSE, cluster_cols = FALSE)
heatmap(df_matrix, scale="column")
heatmap(df_matrix, Colv = NA, Rowv = NA, scale="column", col = coul, xlab="variable", ylab="car", main="heatmap")

######################## NJ TREE #######################

#Realizamos un arbol NJ basado en las distancias FST a partir del conjunto de datos de SNPs filtrados por LD para P.lilfordi
#El metodo de NJ  es un método de agrupación para la creación de árboles fenéticos (fenogramas), creado por Naruya Saitou y Masatoshi Nei en 1987.Se utiliza para árboles de secuencias de ADN o de proteína, para lo cual, el algoritmo requiere del conocimiento de la distancia que existe entre cada par de taxones (por ejemplo, especies o secuencias) para formar el árbol
#En este caso para realizar nuestro NJ usaremos nuestros datos de FST (Diferenciacion genetica)

#Nos ubicamos en la carpeta local desde donde vamos a trabajar, hay que asegúrese de que está en la carpeta correcta con los archivos descargados disponibles.
setwd("/Users/katherin_otalora/Documents/R_analysis/GATK_2022/allSamples/allSamples_0.90")

#Definimos la varible df con una tabla .txt de FST previamente realizada en excel (La misma que se utilizo con el HeatMap), una matriz espejo.

df <- read.table("FST_B10000_allSamples_GATK_22_05_23.txt",header = TRUE)

#Lo convertimos a matriz con la funcion as.matrix

df_matrix <-as.matrix(df, rownames=TRUE, rownames.value=NULL)
df_matrix

#Corremos el analisis de NJ Tree con la funcion nj y usando nuestra matriz (df_matrix) y lo asignamos a una variable nombrada NJ_Plilfordi luego ploteamos con colores y titulo
NJ_Plilfordi <-nj (df_matrix)
NJ_Plilfordi
plot(NJ_Plilfordi, "u", cex = 0.9)
title("NJ Tree Podarcis lilfordi")
add.scale.bar(x, length = NULL, ask = FALSE, lcol = "black", col= "black")
tiplabels () #Optional

##################################Variance explained and BIC K-MEANS #############################

#La DAPC requiere en sí misma la definición de grupos previos. Sin embargo, los grupos son a menudo desconocidos o inciertos, y es necesario identificar los grupos genéticos antes de describirlos.
#Esto puede lograrse mediante k-means, un algoritmo de agrupación que encuentra un número determinado (digamos, k) de grupos que maximiza la variación entre grupos, B(X).
#Para identificar el número óptimo de grupos, k-means se ejecuta secuencialmente con valores crecientes de k, y las diferentes soluciones de agrupación se comparan utilizando el criterio de información bayesiano (BIC).
#Lo ideal es que la solución de agrupación óptima corresponda al BIC más bajo. En la práctica, el "mejor" BIC suele estar indicado por un codo en la curva de valores BIC en función de k.

#gnd.podarcislilfordi es nuestro conjunto de datos de SNPs. Utilizamos find.clusters para identificar los grupos. Especificamos que queremos evaluar hasta k = 16 grupos o 30 grupos (max.n.clust=10), ya que tenemos 16 poblaciones con allSample
#El max.n.clust dependera del numero de localidades o islotes en el caso de P.lilfordi que usted tenga y quiera evaluar, asi como de la posible estructuracion que considera que tiene su data set
#El max.n.clust se puede utilizar el doble de tamano de sus localidades ya que puede que cada localidad este estrcuturada en al menos dos grupos geneticos

grp <- find.clusters(gnd.podarcislilfordi, max.n.clust=30) # fincluster transforma los datos usando PCA, pidiendo al usuario que especifique el número de PC retenidas de forma interactiva a menos que se proporcione el argumento n.pca y ejecuta el algoritmo k-means (función kmeans del paquete stats) con valores crecientes de k, a menos que se proporcione el argumento n.clust, y calcula las estadísticas de resumen asociadas (por defecto, BIC).
#La función muestra un gráfico de la varianza acumulada explicada por los valores propios del PCA. Aparte del tiempo de cálculo, no hay ninguna razón para mantener un número pequeño de componentes; aquí, mantenemos toda la información, especificando que se conserven 200 PC (en realidad hay menos PC -alrededor de 110-, por lo que se mantienen todos).
#Choose the number PCs to retain (>= 1):
250 # Se escoje ese valor porque ya se sabe que tiene menos que 150 PCs aprox por que la pendiente llega a balance en ese numero en la Variance explained
#Luego de este paso smuestra un gráfico de BIC que te dara el numero de clusters a retener y correspondera al valor mas bajo, donde la curva inicia a crecer en este caso 3
#Choose the number of clusters (>=2):
6 # Se escoje ese valor porque es donde cambia drasticamente y comienza a subir la curva del BIC

# RESULTADO BIC 6= 1112.488 

#Este resultado nos sirve para definir cuantos grupos geneticos esperamos tener tanto en el PCA como en el DAPC

##################### CROSS VALIDATION #################
#A continuación, dedicaremos un tiempo a establecer cuál es el número apropiado de componentes principales (PC) para el análisis.
#Es importante elegir cuidadosamente la cantidad correcta de PC para incluir la mayoría de las fuentes de variación explicadas por una cantidad adecuada de PC retenidas.
#Una forma de asegurarse de haber seleccionado la cantidad correcta de PC es realizar un CROSS VALIDATION.
#Este es un procedimiento en el que omite un cierto porcentaje de sus datos, ejecuta DAPC y luego ve si los datos que se omitieron se colocan correctamente en la población.
#Realice el cross validation para encontrar el número óptimo de PCs a retener en el DAPC, vamos a correrlo con 1000 replicas

grp <- pop(gnd.podarcislilfordi)
set.seed(999)
x = tab(gnd.podarcislilfordi, NA.method = "mean")
crossval_1 = xvalDapc(x, gnd.podarcislilfordi$pop, result = "groupMean", n.pca.max = 300, training.set = 0.9, center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 1000, xval.plot = TRUE)

crossval_1

#Aqui vamos a ver los principales resultados para poder seleccionar cuantos PCs debemos usar como numero optimo en el analisis de DAPC
crossval_1[2:6]

#En los resultados nos quedamos con la seccion de:
#$`Number of PCs Achieving Lowest MSE`
#$`Number of PCs Achieving Lowest MSE` #Nos quedamos con 60 para hacer el DAPC
#[1] "60"

#Y tambien con la seccion de:

#$n.pca: 60 first PCs of PCA used
#$n.da: 15 discriminant functions saved
#$var (proportion of conserved variance): 0.632 ###### Este es el porcentaje de varianza que explican los 60 PRIMEROS PCs ###########

######################################### PCA ##################################
#En estadística, el análisis de componentes principales (PCA) es una técnica utilizada para describir un conjunto de datos en términos de nuevas variables («componentes») no correlacionadas.
#Los componentes se ordenan por la cantidad de varianza original que describen, por lo que la técnica es útil para reducir la dimensionalidad de un conjunto de datos.
#Un análisis de componentes principales (PCA) convierte los datos SNPs observados en un conjunto de valores de variables no correlacionadas linealmente llamadas componentes principales que resumen la variación entre las muestras.

#Podemos realizar un PCA en nuestro objeto genlight utilizando la función glPCA y definiendo nf que indica el número de componentes principales que deben retenerse
#Es importante establecer nf = NULL cuando se analizan los datos porque la cantidad de componentes principales retenidos tiene un gran efecto en el resultado de los datos.
#El método estadístico cross validation  ayuda para elegir nf (PCA) y n.pca (DAPC) (Ver script mas arriba para calcular el parametro)
#El nf depende del PCA eigenvalues por eso mejor primero correrlo ver los resultados y volver a correrlo con el nf, tambien esta la opcion interactiva de nf=NULL y luego de ver los Eigenvalues asignar el nf
#El nf se selecciona tratando de escoger el numero de componente principales que explican acumulativamente entre el 60% al 80%de la varianza de los datos.
#Por lo general los 3 primeros componentes principales suelen acumular esa varianza pero en SNPs varia mucho y pueden ser mas de 3


podarcislilfordi.pca <- glPca(gl.podarcislilfordi, nf =60)  #nf 60 por el crossvalidation
#Aqui saldra una grafica y al correlo nos dice si los primero componentes retienen la sumatoria del 80% varianza para poder seleccionarlos
#Como hicimos el analisis de Crossvalidation podemos retener 40
#Choose the number of PCs (>=2):
#3 o mas en este caso 60
barplot(100*podarcislilfordi.pca$eig/sum(podarcislilfordi.pca$eig),col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

#Con esta opcion podemos ver los valores de la varianza que explica cada uno de nuestros componentes dentro del EigenValues
PC1_percetnaje <-  podarcislilfordi.pca.scores %>% apply(MARGIN = 2, FUN = sum)
eigen_percent <- round(((podarcislilfordi.pca$eig)/sum(podarcislilfordi.pca$eig))*100,2)
eigen_percent

[1] 14.53  5.17  3.20  3.07  2.53  1.85  1.66  1.46  1.31  1.20  0.97  0.91
[13]  0.84  0.79  0.78  0.76  0.65  0.64  0.60  0.60  0.58  0.58  0.56  0.56
[25]  0.55  0.55  0.54  0.54  0.53  0.52  0.52  0.52  0.51  0.51  0.50  0.50
[37]  0.50  0.49  0.49  0.49  0.49  0.48  0.48  0.47  0.47  0.47  0.47  0.46
[49]  0.46  0.46  0.46  0.45  0.45  0.45  0.44  0.44  0.44  0.44  0.43  0.43
[61]  0.43  0.43  0.43  0.42  0.41  0.41  0.41  0.41  0.41  0.40  0.40  0.40
[73]  0.40  0.40  0.39  0.39  0.39  0.39  0.38  0.38  0.38  0.38  0.37  0.37
[85]  0.37  0.37  0.36  0.36  0.36  0.36  0.36  0.35  0.35  0.35  0.35  0.35
[97]  0.34  0.34  0.34  0.34  0.33  0.33  0.33  0.33  0.33  0.33  0.32  0.32
[109]  0.32  0.32  0.32  0.31  0.31  0.31  0.31  0.31  0.31  0.30  0.30  0.30
[121]  0.30  0.29  0.29  0.29  0.29  0.29  0.28  0.28  0.28  0.28  0.28  0.28
[133]  0.27  0.27  0.27  0.27  0.27  0.26  0.26  0.26  0.26  0.26  0.25  0.25
[145]  0.25  0.25  0.25  0.24  0.24  0.24  0.24  0.24  0.23  0.23  0.23  0.23
[157]  0.23  0.23  0.22  0.22  0.22  0.21  0.21  0.21  0.21  0.20  0.20  0.20
[169]  0.20  0.19  0.19  0.19  0.18  0.18  0.17  0.17  0.16  0.16  0.15  0.15
[181]  0.15  0.15  0.14  0.14  0.13  0.12  0.11  0.06  0.04  0.02

#Para ver los resultados del PCA podemos utilizar el paquete ggplot2. Tenemos que convertir los datos que contiene los componentes principales (podarcislilfordi.pca$scores) en el nuevo objeto podarcislilfordi.pca.scores.
#Además, añadiremos los valores de la población como una nueva columna en nuestro objeto podarcislilfordi.pca.scores, para poder colorear las muestras por población.

podarcislilfordi.pca.scores <- as.data.frame(podarcislilfordi.pca$scores)
podarcislilfordi.pca.scores$pop <- pop(gl.podarcislilfordi)

p <- ggplot(podarcislilfordi.pca.scores, aes(PC1, PC2, colour=pop))
p <- p + geom_point(size=1)
p <- p + labs(x = "PC1 (14.53%)", y = "PC2 (5.17%)") #En esta seccion hay que cambiarle los valores de cada PC despendiendo del resultado del eigen_percent (ver script mas arriba)
p <- p + stat_ellipse(level =0,99, size = 1)
p <- p + coord_equal() + theme_light()
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + theme_bw()

p

################ DAPC and COMPOPLOT  #####################
#DAPC es un enfoque estadístico multivariante que utiliza poblaciones definidas a priori para maximizar la varianza entre las poblaciones de la muestra, dividiéndola en componentes entre poblaciones y dentro de ellas.
#Así, la DAPC maximiza la discriminación entre grupos.
#DAPC requiere un objeto genlight con poblaciones definidas a priori.
#Nosotros ya tenemos este objeto genlight desde los pasos anteriores. Normalmente, utilizamos el número de componentes principales y los ejes discriminantes que maximizan la varianza entre las poblaciones; pero nuestro objetivo aquí es calcular las asignaciones de población basadas en los resultados del PCA.
#Utilizaremos los mismos parámetros que en el PCA para que los resultados sean comparables entre ambos métodos. Estos parámetros (n.pca=20 y n.da=9) se utilizarán para reconstruir el DAPC, obtener la asignación de las muestras a cada población y sugerir la mezcla entre localidades o islotes.
#Es importante establecer n.pca = NULL cuando se analizan los datos y no se ha hecho el crossvalidation porque la cantidad de componentes principales retenidos tiene un gran efecto en el resultado de los datos.
#El método estadístico cross validation  ayuda para elegir n.pca (Ver script mas arriba para calcular el parametro)
#El n.pca= sera igual al resultado obtenido por Crossvalidation
# El n.da= es el número de ejes retenidos en el Análisis Discriminante (DA) tambien se selecciona con el Crossvalidation

pnw.dapc <- dapc(gl.podarcislilfordi, n.pca = 60, n.da =15)
nb.cols <- 16
cols <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
scatter(pnw.dapc, col =cols, cex =2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.8, ratio.da=.18, ratio.pca=.18)
#title(ylab="LD2", line = 2) #opcional
#title(xlab="LD1", line = 1) #opcional

# HACER EL GRAFICO Con puntos transparentes
scatter(pnw.dapc, col =cols, cex =1.8, solid=.1, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.8, ratio.da=.18, ratio.pca=.18)

#Otras opciones que no se dejaron para el paper
#scatter(pnw.dapc, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=cols, solid=.2, cex=3, clab=0, leg=TRUE)
#scatter(pnw.dapc,1,1, col=cols, bg="white",scree.da=TRUE, legend=TRUE, solid=.4)
#scatter(pnw.dapc,2,1, col =cols, cex =1.8, solid=.1, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.8, ratio.da=.18, ratio.pca=.18)

#Con el siguiente script podemos analizar el porcentaje de varianza genética que explica cada eje
percent = pnw.dapc$eig/sum(pnw.dapc$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,100),
        names.arg = round(percent, 1))
#percent
#[1] 34.9955015 20.2446837  9.8811375  8.9086925  6.4037632  3.7871427
#[7]  3.5304288  3.0607997  2.4281206  2.3059897  1.7749455  1.0370859
#[13]  0.9557512  0.3964812  0.2894763


#Hacemos el COMPOPLOT gráfico de barras apiladas compuestas
#El objeto DAPC que creamos incluye la probabilidad de pertenencia a la población de cada muestra para cada una de las poblaciones predeterminadas. Para visualizar las asignaciones posteriores de cada muestra, utilizamos un gráfico de barras apiladas compuestas (compoplot).
#Un complot ilustra la probabilidad de pertenencia a la población en el eje y.
#Cada muestra es un "BIN" en el eje x, y la probabilidad asignada de pertenencia a la población se muestra como un gráfico de barras apiladas con conglomerados o poblaciones en color

compoplot(pnw.dapc,col = cols, posi = 'bottom')

#Estos gráficos son difíciles de interpretar, por lo que separaremos las muestras por población, para hacerlo mas sencillo.
#Se puede utilizar ggplot2 para reconstruir estos gráficos, pero tenemos que convertir los datos en un objeto de ggplot2.
#Extraeremos las asignaciones de pertenencia a la población calculadas por DAPC (pnw.dapc$posterior) en un nuevo marco de datos (dapc.results), incluiremos la asignación original de la población como una nueva columna en los datos (dapc.results$pop)
#Y añadiremos una columna que incluya los nombres de las muestras (dapc.results$indNames).

dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.podarcislilfordi)
dapc.results$indNames <- rownames(dapc.results)

#ggplot2 tiene requisitos específicos para la estructura del formato de datos, ya que requiere cada observación en filas, y todos los diferentes valores de estas observaciones en columnas. Para transformar los datos utilizamos pivot_longer del paquete tidyr.

dapc.results <- melt(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

#Ahora hemos reorganizado los datos en el formato requerido, donde cada observación de probabilidad de pertenencia para una población dada es una fila con el nombre de la muestra, la población original y la población asignada como columnas.

head(dapc.results, n = 6)

#A continuación, cambiamos el nombre de las columnas a términos más familiares:

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

#ggplot2 graficara los datos dapc.results que hemos reorganizado utilizando pivot_longer, utilizando las muestras en el eje X y las probabilidades de pertenencia en el eje Y. El color de relleno indicará las asignaciones poblacionales originales. Cada "facet" representa la asignación original de la población para cada muestra:

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity')
p <- p + scale_fill_manual(values = cols)
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
p

#Otra opcion para ver la dispersión de los datos
#Hay que tener en cuenta que la dispersión también puede representar una única función discriminante, lo que resulta especialmente útil cuando sólo se ha conservado una de ellas (por ejemplo, en el caso k = 2).
#Esto se consigue trazando las densidades de los individuos en una función discriminante determinada con colores diferentes para los distintos grupos:
scatter(pnw.dapc,1,1, col=cols, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.6)
scatter(pnw.dapc,2,2, col=cols, bg="white",
                scree.da=FALSE, legend=TRUE, solid=.6)
#resultados

pnw.dapc
#################################################
# Discriminant Analysis of Principal Components #
#################################################
class: dapc
$call: dapc.genlight(x = gl.podarcislilfordi, n.pca = 60, n.da = 15)

$n.pca: 60 first PCs of PCA used
$n.da: 15 discriminant functions saved
$var (proportion of conserved variance): 0.632

$eig (eigenvalues): 5201 3009 1468 1324 951.7 ...

vector    length content                   
1 $eig      15     eigenvalues               
2 $grp      191    prior group assignment    
3 $prior    16     prior group probabilities 
4 $assign   191    posterior group assignment
5 $pca.cent 2876   centring vector of PCA    
6 $pca.norm 2876   scaling vector of PCA     
7 $pca.eig  190    eigenvalues of PCA        

data.frame    nrow ncol content                                          
1 $tab          191  60   retained PCs of PCA                              
2 $means        16   60   group means                                      
3 $loadings     60   15   loadings of variables                            
4 $ind.coord    191  15   coordinates of individuals (principal components)
5 $grp.coord    16   15   coordinates of groups                            
6 $posterior    191  16   posterior membership probabilities               
7 $pca.loadings 2876 60   PCA loadings of original variables               
8 $var.contr    2876 15   contribution of original variables  


##################### SUBSET  by population #################################

popNames(gl.podarcislilfordi)

Menorca <- popsub(gl.podarcislilfordi, sublist=c(1:2, 5, 13:16))
popNames(Menorca)

Mallorca <-popsub(gl.podarcislilfordi, sublist=c(3:4, 6:12))
popNames(Mallorca)


#### MENORCA ####

#### CROSSVALIDATION MENORCA #####
grp <- pop(Menorca)
set.seed(999)
x = tab(Menorca, NA.method = "mean")
crossval_Menorca = xvalDapc(x, Menorca$pop, result = "groupMean", n.pca.max = 300, training.set = 0.9, center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 1000, xval.plot = TRUE)

crossval_Menorca

#Aqui vamos a ver los principales resultados para poder seleccionar cuantos PCs debemos usar como numero optimo en el analisis de DAPC
crossval_Menorca[2:6]
#$`Number of PCs Achieving Lowest MSE`
#[1] "40"

#### FINAL DAPC with 40 PCs MENORCA ####
#now let's redo the dapc analyses with the correct n.pca to retained
dapcfinal_Menorca <- dapc(Menorca, var.contrib = TRUE, scale = FALSE, n.pca = 40, n.da = nPop(Menorca) - 1) 
nb.cols <- nPop(Menorca) - 1
cols <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
scatter(dapcfinal_Menorca, col =cols, cex =2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.8, ratio.da=.18, ratio.pca=.18)

#Hacemos el COMPOPLOT gráfico de barras apiladas compuestas
#El objeto DAPC que creamos incluye la probabilidad de pertenencia a la población de cada muestra para cada una de las poblaciones predeterminadas. Para visualizar las asignaciones posteriores de cada muestra, utilizamos un gráfico de barras apiladas compuestas (compoplot).
#Un complot ilustra la probabilidad de pertenencia a la población en el eje y.
#Cada muestra es un "BIN" en el eje x, y la probabilidad asignada de pertenencia a la población se muestra como un gráfico de barras apiladas con conglomerados o poblaciones en color

compoplot(dapcfinal_Menorca,col = cols, posi = 'bottom')

#Estos gráficos son difíciles de interpretar, por lo que separaremos las muestras por población, para hacerlo mas sencillo.
#Se puede utilizar ggplot2 para reconstruir estos gráficos, pero tenemos que convertir los datos en un objeto de ggplot2.
#Extraeremos las asignaciones de pertenencia a la población calculadas por DAPC (pnw.dapc$posterior) en un nuevo marco de datos (dapc.results), incluiremos la asignación original de la población como una nueva columna en los datos (dapc.results$pop)
#Y añadiremos una columna que incluya los nombres de las muestras (dapc.results$indNames).

dapc.results <- as.data.frame(dapcfinal_Menorca$posterior)
dapc.results$pop <- pop(Menorca)
dapc.results$indNames <- rownames(dapc.results)

#ggplot2 tiene requisitos específicos para la estructura del formato de datos, ya que requiere cada observación en filas, y todos los diferentes valores de estas observaciones en columnas. Para transformar los datos utilizamos pivot_longer del paquete tidyr.

dapc.results <- melt(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

#Ahora hemos reorganizado los datos en el formato requerido, donde cada observación de probabilidad de pertenencia para una población dada es una fila con el nombre de la muestra, la población original y la población asignada como columnas.

head(dapc.results, n = 6)

#A continuación, cambiamos el nombre de las columnas a términos más familiares:

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

#ggplot2 graficara los datos dapc.results que hemos reorganizado utilizando pivot_longer, utilizando las muestras en el eje X y las probabilidades de pertenencia en el eje Y. El color de relleno indicará las asignaciones poblacionales originales. Cada "facet" representa la asignación original de la población para cada muestra:

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity')
p <- p + scale_fill_manual(values = cols)
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
p

#### MALLORCA ####

#### CROSSVALIDATION MENORCA #####
grp <- pop(Mallorca)
set.seed(999)
x = tab(Mallorca, NA.method = "mean")
crossval_Mallorca = xvalDapc(x, Mallorca$pop, result = "groupMean", n.pca.max = 300, training.set = 0.9, center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 1000, xval.plot = TRUE)

crossval_Mallorca

#Aqui vamos a ver los principales resultados para poder seleccionar cuantos PCs debemos usar como numero optimo en el analisis de DAPC
crossval_Mallorca[2:6]
#$`Number of PCs Achieving Lowest MSE`
#[1] "30"

#### FINAL DAPC with 30 PCs MALLORCA ####
#now let's redo the dapc analyses with the correct n.pca to retained
dapcfinal_Mallorca <- dapc(Mallorca, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(Mallorca) - 1) 
nb.cols <- nPop(Mallorca) - 1
cols <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
scatter(dapcfinal_Mallorca, col =cols, cex =2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topright", cleg = 0.8, ratio.da=.18, ratio.pca=.18)

#Hacemos el COMPOPLOT gráfico de barras apiladas compuestas
#El objeto DAPC que creamos incluye la probabilidad de pertenencia a la población de cada muestra para cada una de las poblaciones predeterminadas. Para visualizar las asignaciones posteriores de cada muestra, utilizamos un gráfico de barras apiladas compuestas (compoplot).
#Un complot ilustra la probabilidad de pertenencia a la población en el eje y.
#Cada muestra es un "BIN" en el eje x, y la probabilidad asignada de pertenencia a la población se muestra como un gráfico de barras apiladas con conglomerados o poblaciones en color

compoplot(dapcfinal_Menorca,col = cols, posi = 'bottom')

#Estos gráficos son difíciles de interpretar, por lo que separaremos las muestras por población, para hacerlo mas sencillo.
#Se puede utilizar ggplot2 para reconstruir estos gráficos, pero tenemos que convertir los datos en un objeto de ggplot2.
#Extraeremos las asignaciones de pertenencia a la población calculadas por DAPC (pnw.dapc$posterior) en un nuevo marco de datos (dapc.results), incluiremos la asignación original de la población como una nueva columna en los datos (dapc.results$pop)
#Y añadiremos una columna que incluya los nombres de las muestras (dapc.results$indNames).

dapc.results <- as.data.frame(dapcfinal_Menorca$posterior)
dapc.results$pop <- pop(Menorca)
dapc.results$indNames <- rownames(dapc.results)

#ggplot2 tiene requisitos específicos para la estructura del formato de datos, ya que requiere cada observación en filas, y todos los diferentes valores de estas observaciones en columnas. Para transformar los datos utilizamos pivot_longer del paquete tidyr.

dapc.results <- melt(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

#Ahora hemos reorganizado los datos en el formato requerido, donde cada observación de probabilidad de pertenencia para una población dada es una fila con el nombre de la muestra, la población original y la población asignada como columnas.

head(dapc.results, n = 6)

#A continuación, cambiamos el nombre de las columnas a términos más familiares:

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

#ggplot2 graficara los datos dapc.results que hemos reorganizado utilizando pivot_longer, utilizando las muestras en el eje X y las probabilidades de pertenencia en el eje Y. El color de relleno indicará las asignaciones poblacionales originales. Cada "facet" representa la asignación original de la población para cada muestra:

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity')
p <- p + scale_fill_manual(values = cols)
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
p


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


############################# Divmigrate ###########################################

##Crear objeto Genepop## usando el programa PGDspider
#El output de PGDspider es .txt por eso hay que abrir el archivo y ahora guardarlo con la extension .gen
#OJO IMPORTANTE PONER POP cuando inicia una nueva poblacion para delimitarlas, si no saldra el error de:"" Error in combn(npops, 2) : n < m"
#Para los datos de extSample tuve que organizar los individuos de aire y ubicarlos al final que esta su poblacion, y luego correr el analisis.

#Otra alternativa:
#gnd.podarcislilfordi_genpop <- genind2genpop(gnd.podarcislilfordi)
#gnd.podarcislilfordi_genpop

#INSTALAR Y LLAMAR LAS LIBRERIAS QUE SE NECESITAN
install.packages("DNAprofiles")
library(DNAprofiles)
install.packages("diveRsity")
library(diveRsity)

setwd("/Users/katherin_otalora/Documents/Divmigrate/GENEPOP")
source("/Users/katherin_otalora/Downloads/diveRsity/man/divMigrate.Rd")

#El output de PGDspideres .txt por eso hay que abrir el archivo y ahora guardarlo con la extension .gen
Resuts_migration_plilfordi <- divMigrate(infile ="/Users/katherin_otalora/Documents/Divmigrate/GENEPOP/allSample_0.90_filtredLD_genepop_1_1.gen", boots = 1000, stat = "gst",plot_network = TRUE,  plot_col = "darkblue", para =8)

#Ahora solo mostraremos el grafico con los valores superiores o iguales a 0,5 de tasa de migracion, no es necesario hacer el boots porque siempre la matriz de flujo y el boots da igual ########
Resuts_migration_plilfordi <- divMigrate(infile ="/Users/katherin_otalora/Documents/Divmigrate/GENEPOP/allSample_0.90_filtredLD_genepop_1_1.gen", boots = 1000, stat = "gst",filter_threshold = 0.25, plot_network = TRUE,  plot_col = "darkblue", para =8)
Resuts_migration_plilfordi <- divMigrate(infile ="/Users/katherin_otalora/Documents/Divmigrate/GENEPOP/allSample_0.90_filtredLD_genepop_1_1.gen", boots = 1000, stat = "gst",filter_threshold = 0.5, plot_network = TRUE,  plot_col = "darkblue", para =8)

#Salen dos graticos uno con el Relative migration network y posteriormente el significat relative migration network

##### Por Islas ######

Resuts_migration_plilfordi_island <- divMigrate(infile ="/Users/Katherin_Otalora/Desktop/Katherin_PhD/DivMigrate/raw.plilfordi_run05_03_04_maf_005_hwe_indels_maxmiss_0_02_1SNP-locus_Mallorca_Menorca", outfile = "Plilfordi_Menoca_Mallorca_0.25", stat = "gst",filter_threshold = 0.25 , plot_network = TRUE,  plot_col = "darkblue", para =TRUE)

Resuts_migration_plilfordi_island <- divMigrate(infile ="/Users/Katherin_Otalora/Desktop/Katherin_PhD/DivMigrate/raw.plilfordi_run05_03_04_maf_005_hwe_indels_maxmiss_0_02_1SNP-locus_Mallorca_Menorca", outfile = "Plilfordi_Menoca_Mallorca_0.25", boots = 1000, stat = "gst",filter_threshold = 0.25 , plot_network = TRUE,  plot_col = "darkblue", para =TRUE)

