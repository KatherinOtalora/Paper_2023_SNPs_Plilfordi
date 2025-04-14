############################################################# GEA GENOME ENVIROMENT ASSOCIATION ###########################################################################
#https://popgen.nescent.org/2018-03-27_RDA_GEA.html

# Load packages
# -------------
library(vcfR)
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA
library(adegenet)
library(dplyr)

#Nos ubicamos en la carpeta local desde donde vamos a trabajar, hay que asegúrese de que está en la carpeta correcta con los archivos descargados disponibles.

setwd("/Users/katherin_otalora/Documents/GEA/GEA_extsamples_ouliers_09_06_23")

#Read VCF and transform it to a genotype matrix
vcf <- read.vcfR("outliers_filtred_LOGPO_FST_comunes_bayescan_pcadapt_extsamples_withoutLD_1_100_K=10.recode.vcf", verbose = FALSE )
gt <- extract.gt(vcf, element = c('GT'), as.numeric = TRUE)
head(gt)
#Veremos las dimensiones de este archivo
dim(gt)
gen=t(gt) #trasnpose the matrix
#190  91

# NAs
#Calculamos el porcentaje de datos que faltan. Los análisis multivariantes son muy sensibles a los datos que faltan y hay que rellenarlos.
sum(is.na(gt))
#1069

#Los valores que faltan se rellenan utilizando la media global de la variación de la frecuencia alélica.
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

#Luego volvemos a calcular el porcentaje de datos faltantes y debería dar cero.
sum(is.na(gen.imp)) # No NAs
#0

#Descargamos los datos ambientales en R.
#La tabla de las variables se arreglo en excel
#Read the environmental variables
env  <- read.csv (file='Prueba_variables_ambientales_2.csv',sep=";",)

str(env) # Look at the structure of the data frame

#Ahora creamos variables cambiandolas a caracteres o factores dependiendo de cada variable
env$AccessID <- as.character(env$AccessID) #Hacer deL AccessID sean caracteres (no factores)
env$Population <- as.character(env$Population) # Haz que las Population de los individuos sean caracteres (no factores)

str(env)

# Confirmar que los genotipos y los datos ambientales están en el mismo orden

identical(rownames(gen.imp), env[,1])

#El RDA es un método basado en la regresión, por lo que puede estar sujeto a problemas cuando se utilizan predictores altamente correlacionados (Dormann et al., 2013).
#En general, la "regla empírica" |r| > 0,7 es una buena directriz para eliminar los predictores correlacionados. También comprobaremos la multicolinealidad utilizando los Factores de Inflación de la Varianza (VIF), más adelante.

#Evaluamos la correlación entre nuestras variables
#Cambiar el numero despues del : segun el numero de variables que tengamos en este caso 8
pairs.panels(env[,3:5], scale=T)

#Superficie y altura están correlacionadas a 0,75, y Capacidad biotica y altura 0.77.
#Todas tienen valores muy altos de correlacion
#Asumimos el error y no eliminamos ninguna variable? o eliminamos Capacidad biotica?
#NO ELIMINAMOS NINGUNA
#pred <- env[,3:5]
#str(pred)
#pairs.panels(pred, scale=T)

#En caso de querer eliminar alguna usar el siguiente script
#Eliminamos Altitude
pred_1 <- subset(env, select=-c(Maximum_Altitude))
str(pred_1)
pred_1 <- pred_1[,3:4]
str(pred_1)
pairs.panels(pred_1, scale=T)

#En caso de querer eliminar alguna usar el siguiente script
#Extraemos altitude? o capacidad biotica, y maxima altitud, para quedarnos con el resto de las variables
#pred_1 <- subset(env, select=-c(Maximum_Altitude))
#pred_1 <- subset(env, select=-c(Surface))
#pred <- subset(pred_1, select=-c(Biotic_capacity))

#Extraemos altitude, para quedarnos con el resto de las variables
#pred_1 <- subset(env, select=-c(Maximum_Altitude))
#str(pred_1)

#pred_1 <- subset(env)
#str(pred_1)

#Es importante ver la distribución de los niveles de un factor con la siguiente opción #table (pred$surface), ya que en algunos casos no se distribuye equitativamente es decir es sesgada y no hay homogenizacion de la información.
#Veamos este conjunto reducido de predictores:
#pred <- pred_1[,c(3,7)]

#pred <- pred_1[,c(3,4,8)]
#str(pred)


##################################### Ejecutar el análisis de redundancia (RDA) ##################################################

#Este paso puede llevar un tiempo dependiendo de su conjunto de datos.

plilfordi.rda <- rda(gen.imp ~ ., data=pred_1, scale=T)

summary(plilfordi.rda)

#Inertia is variance in the output 
#El R^2 informa sobre el porcentaje de variación genómica que puede ser explicado por uno de los predictores. Ejemplo: cuando encontramos un R^2 ajustado de 0,0002, significa que nuestra ordenación restringida explica alrededor del 0,02% de la variación. Este bajo poder explicativo no es sorprendente dado que esperamos que la mayoría de los SNPs en nuestro conjunto de datos no muestren una relación con los predictores ambientales (por ejemplo, la mayoría de los SNPs serán neutrales, si son outliers esperamos mas asociacion).

#En nuestro caso esperaríamos un R^2 alto para poder ver explicación con nuestros outliers que no serian SNPs neutrales.

RsquareAdj(plilfordi.rda)

#$r.squared
#[1] 0.2586093

#$adj.r.squared
#[1] 0.2417595

#Nuestro R^2 dio 0.34 es decir un 34% de la variación. Este es un buen poder explicativo de la variación y lo esperabamos asi, ya que nuestros datos son outliers y deberían mostrar una relación con los predicadores ambientales.

#Observamos que los valores propios de los ejes restringidos reflejan la varianza explicada por cada eje canónico.

summary(eigenvals(plilfordi.rda, model = "constrained"))


#Podemos visualizar esta información utilizando un gráfico de escala de los valores propios canónicos llamando a screeplot:

screeplot(plilfordi.rda)


#Aquí, podemos ver que el primer eje explica la mayor parte de la varianza. El screeplot proporciona una manera informal (y rápida) de determinar cuántos ejes restringidos incluir cuando buscamos SNPs candidatos. Podríamos empezar investigando los ejes RDA que explican la mayor parte de la varianza (excluyendo los que están después del punto de "caída" en el screeplot).

#Ahora vamos a comprobar la significación de nuestro modelo RDA mediante pruebas formales. Podemos evaluar tanto el modelo completo como cada eje restringido utilizando los estadísticos F (Legendre et al, 2010). La hipótesis nula es que no existe ninguna relación lineal entre los datos del SNP y los predictores ambientales. Consulte anova.cca para obtener más detalles y opciones.

# default is permutation=999
signif.full <- anova.cca(plilfordi.rda, parallel=getOption("mc.cores"), permutations = how(nperm=1000))
signif.full

#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 1000

#Model: rda(formula = gen.imp ~ Surface + Biotic_capacity, data = pred_1, scale = T)
#Df Variance      F   Pr(>F)    
#Model     2   49.136 15.348 0.000999 ***
#  Residual 88  140.864                    
#---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#El modelo completo es significativo, pero eso no nos dice mucho. 
#Podemos comprobar la significación de cada eje restringido utilizando el siguiente código.
#Este análisis requiere mucho tiempo (hasta unas horas dependiendo del conjunto de datos). 

signif.axis <- anova.cca(plilfordi.rda, by="axis", parallel=getOption("mc.cores"), permutations = how(nperm=1000))
signif.axis

#Permutation test for rda under reduced model
#Forward tests for axes
#Permutation: free
#Number of permutations: 1000

#Model: rda(formula = gen.imp ~ Surface + Biotic_capacity, data = pred_1, scale = T)
#Df Variance       F   Pr(>F)    
#RDA1      1   47.470 29.6549 0.000999 ***
#  RDA2      1    1.666  1.0409 0.334665    
#Residual 88  140.864                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Encontramos que el primer ejes restringido es significativos (p = 0.000999); el eje restringido 2 tiene un valor p de 0,3346, mientras que los eje 3 tienen valores p =0.7642. Esto se corresponde con nuestra evaluación del diagrama de escala, más arriba.
#En nuestros resultados solo el RDA1 fue significativo

#Por último, vegan dispone de una función sencilla para comprobar los factores de inflación de la varianza de las variables predictoras utilizadas en el modelo:

vif.cca(plilfordi.rda)

#Surface        Biotic_capacity 
#2.496557           2.496557 

#Todos los valores son inferiores a 10, y la mayoría son inferiores a 5, lo que indica que la multicolinealidad entre estos predictores no debería ser un problema para el modelo. Podríamos eliminar una de las variables X si nos preocuparan estos valores VIF más altos (Zuur et al., 2010).

############################### PLOT RDA #############################################################

# default is axes 1 and 2

plot(plilfordi.rda, scaling=3)

plot(plilfordi.rda, choices = c(1, 3), scaling=3)  # axes 1 and 3

#En el gráfico, los SNPs están en rojo (en el centro de cada gráfico), y los individuos son los círculos negros. 
#Los vectores azules son los predictores ambientales. La disposición relativa de estos elementos en el espacio de ordenación refleja su relación con los ejes de ordenación, que son combinaciones lineales de las variables predictoras.

#Hagamos algunos gráficos más informativos. Colorearemos los puntos individuales en función de su islote, que podemos encontrar en el conjunto de datos env

levels(env$Population) <- c("REI","AIRE","PORROS","FORADADA","ESCLTS","DRAGONERA","CABRERA","COLOM", "COLOMER")
eco <- env$Population

#Seleccionamos una paleta de colores para cada islote

bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#1adce3","#e31a1c","#cce31a","#781ae3","#d6e31a") # 8 nice colors for our ecotypes

#Esta vez, vamos a configurar los gráficos y añadir cada componente por separado:

plot(plilfordi.rda, type="n", scaling=3)
points(plilfordi.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(plilfordi.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the wolves
text(plilfordi.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#Extraemos la información de los individuos y la guardamos para su publicación.

rda_indv <- as.data.frame(scores(plilfordi.rda, display=c("sites")))
write.table(rda_indv, "rda_plilfordi_91ind_extsamples_outliers_variables_2.txt", quote=FALSE, row.names=TRUE)


########################################### Identificar los SNPs candidatos implicados en la adaptación local ########################################

#Usaremos las cargas de los SNPs en el espacio de ordenación para determinar qué SNPs son candidatos a la adaptación local.
#Las cargas de los SNPs se almacenan como especies en el objeto RDA. Extraeremos las cargas de los SNPs de los tres ejes significativos restringidos

load.rda <- scores(plilfordi.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
#load.rda <- scores(plilfordi.rda, choices=c(1), display="species")  # Solo el RDA1 Species scores for the first three constrained axes

#Si observamos los histogramas de las cargas en cada eje de la RDA, podemos ver sus distribuciones (no son normales ya que estamos trabajando con outliers). Los SNPs que se cargan en el centro de la distribución no muestran una relación con los predictores ambientales; los que se cargan en las colas sí, y es más probable que estén bajo selección en función de esos predictores (o de algún otro predictor correlacionado con ellos).
#Vemos el histograma y SNPs de las colas tienen relacion con las variables ambientales
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")

#Se ha escrito una función sencilla para identificar los SNPs que se cargan en las colas de estas distribuciones. Empezaremos con un límite de 3 desviaciones estándar (valor p de dos colas = 0,0027). Como con todos los límites, esto puede modificarse para reflejar los objetivos del análisis y nuestra tolerancia a los verdaderos positivos frente a los falsos positivos. Por ejemplo, si se necesita ser muy conservador y sólo identificar aquellos loci bajo una selección muy fuerte (es decir, minimizar las tasas de falsos positivos), se podría aumentar el número de desviaciones estándar a 3,5 (valor p de dos colas = 0,0005). Esto también aumentaría la tasa de falsos negativos. Si le preocupan menos los falsos positivos y más la identificación de tantos loci candidatos potenciales como sea posible (incluidos los que pueden estar bajo una selección más débil), podría elegir un límite de desviación estándar de 2,5 (valor p de dos colas = 0,012).

#Aquí defino la función como valores atípicos, donde x es el vector de cargas y z es el número de desviaciones estándar a utilizar:

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#Ahora vamos a aplicarlo a cada uno de los ejes restringidos significativos (0.0093 diferencia entre c/u):
#limite de 3,5 desviaciones estandar (valor p de dos colas=0,0005)
#límite de 3 desviaciones estándar (valor p de dos colas = 0,0027)
#límite de 2,5  desviaciones estándar (valor p de dos colas = 0,012)
#Nosotros decidimos ponerle límite de desviación estándar de 2 (valor p de dos colas= 0.0213)
#límite de 1,5  desviaciones estándar (valor p de dos colas = 0,0306)
#límite de 1  desviaciones estándar (valor p de dos colas = 0,0399)
#límite de 0.5  desviaciones estándar (valor p de dos colas = 0,0492)


cand1 <- outliers(load.rda[,1],0.5) # 46
length(cand1)
cand2 <- outliers(load.rda[,2],0.5) #154
length(cand2)
cand3 <- outliers(load.rda[,3],0.5) #78
length(cand3)

#ncand <- length(cand1)
#ncand

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand
#278

#A continuación, organizaremos nuestros resultados haciendo un marco de datos con el eje, el nombre del SNP, la carga y la correlación con cada predictor:

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
#colnames(cand1) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
#cand <- rbind(cand1)

write.table(cand, "snps_outliers_plilfordi_extsamples_variables_2.txt", quote=FALSE, row.names=TRUE)


#Añadamos las correlaciones de cada SNP candidato con los predictores ambientales:

#foo <- matrix(nrow=(ncand), ncol=3)  # 3 columns for 3 predictors
#colnames(foo) <- c("Surface","Maximum_Altitude", "Biotic_capacity")

foo <- matrix(nrow=(ncand), ncol=2)  # 2 columns for 2 predictors
colnames(foo) <- c("Surface", "Biotic_capacity")


#for (i in 1:length(cand$snp)) {
# nam <- cand[i,2]
#snp.gen <- gen.imp[,nam]
#snp.mat=as.matrix(snp.gen)
#foo[i,] <- apply(pred, 2, function (x) cor(x,snp.mat))
#}

#Este script para 2 predictores con pred_1
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred_1,2,function(x) cor(x,snp.gen))
}


cand <- cbind.data.frame(cand,foo)  
head(cand)

write.table(cand, "correlation_snps_outliers_plilfordi_extsamples_variables_2.txt", quote=FALSE, row.names=TRUE)

###Investigate the candidates####

length(cand$snp[duplicated(cand$snp)])  # 107 duplicate detections

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) #  0 duplicates on axis 1
## 
##  0 
## 46

table(foo[foo[,1]==2,2]) #  0 duplicates on axis 2

#0   1
#125  29

table(foo[foo[,1]==3,2]) # 1 duplicates on axis 3

#1 
#78

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
length(cand$snp)
#171

write.table(cand, "NEW_correlation_snps_outliers_plilfordi_extsamples_variables_2.txt", quote=FALSE, row.names=TRUE)


#We’ve now reduced our candidates to 171 unique SNPs.

#Este script cuando se usen las 3 variables

#for (i in 1:length(cand$snp)) {
# bar <- cand[i,]
# cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable
# cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation
#}

#colnames(cand)[7] <- "predictor"
#colnames(cand)[8] <- "correlation"

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,6] <- names(which.max(abs(bar[4:5]))) # gives the variable
  cand[i,7] <- max(abs(bar[4:5]))              # gives the correlation
}


colnames(cand)[6] <- "predictor"
colnames(cand)[7] <- "correlation"

table(cand$predictor) 

snps_finales <- table(cand$predictor) 
snps_finales 

#Surface 
#171 

write.table(cand, "predictor_ambiental_snps_outliers_plilfordi_extsamples_variables_2.txt", quote=FALSE, row.names=TRUE)


###Plot the SNPs###

sel <- cand$snp
env <- cand$predictor
env[env=="Biotic_capacity"] <- '#1f78b4'
env[env=="Surface"] <- '#a6cee3'
env[env=="Maximum_Altitude"] <- '#6a3d9a'
env

# color by predictor:
col.pred <- rownames(plilfordi.rda$CCA$v) # pull the SNP names
col.pred
for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("Chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
col.pred[grep("scaffold",col.pred)] <- '#f1eef6' # non-candidate SNPs

empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3')

# axes 1 & 2
plot(plilfordi.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(plilfordi.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(plilfordi.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(plilfordi.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("Biotic_capacity","Surface","Maximum_Altitude"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)


# axes 1 & 3
plot(plilfordi.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(plilfordi.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(plilfordi.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(plilfordi.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("Biotic_capacity","Surface", "Maximum_Altitude"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)


