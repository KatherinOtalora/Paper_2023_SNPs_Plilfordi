############ MANTEL TEST ##########
######################################### Allsamples (16 islotes) ##########################################################################
#Primero hacemos calculamos distancia euclidean y geodesica
#Importante: Tienes una lista de coordenadas geográficas medidas con latitud y longitud. 
#Estas coordenadas se miden en grados, y una distancia de 1 grado (especialmente la longitud) no equivale a la misma distancia real (en metros) dependiendo de dónde se encuentre en el globo (mucho mayor en el ecuador que en los polos).
#Se calculan distancias euclideanas y geodésicas para ver la diferencia

install.packages ("geodist")
library(geodist)
setwd("/Users/katherin_otalora/Documents/Scripts_R/FINAL_SCRIPT_ALLSAMPLES/MANTEL")
#coord_dec1 <- read.table("Coordenadas_decimales_Plilfordi_30_05_23.txt", header = TRUE)
coord_dec <- read.table("1Number_Coordenadas_decimales_Plilfordi_30_05_23.txt", header = TRUE)
coord_dec

#coord_dec1
#euclidean distances in degrees:
euclidean_distances <- dist(coord_dec,diag=T, upper=T)
euclidean_distances
euclidean_distances_matrix  <-as.matrix(euclidean_distances, rownames=TRUE, rownames.value=NULL)
euclidean_distances_matrix 

#dist(coord_dec1, diag=T, upper=T)
#geodesic distances: 
geo_distances_matrix <-geodist(coord_dec)
geo_distances_matrix
#SALE UN WARNING PREGUNTAR??

#Se lee la matriz de datos Geneticos
setwd("/Users/katherin_otalora/Documents/Scripts_R/FINAL_SCRIPT_ALLSAMPLES/MANTEL")
#df1 <- read.table("FST_Plilfordi_30_05_23.txt", header = TRUE)
df <- read.table("1FST_Plilfordi_30_05_23.txt", header = TRUE)
df
#df1
genetic_matrix <-as.matrix(df, rownames=TRUE, rownames.value=NULL)
#genetic_matrix <-as.matrix(df1, rownames=TRUE, rownames.value=NULL)
genetic_matrix

#Se realiza el TEST DE MANTEL usando primero la distancia euclideana (Allsamples 16 islotes)
#install.packages ("vegan")
library (vegan)

mantel.p_lilfordi_euc<-mantel(euclidean_distances_matrix,genetic_matrix,method="pearson",permutations=10000)
mantel.p_lilfordi_euc

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#mantel(xdis = euclidean_distances_matrix, ydis = genetic_matrix,      method = "pearson", permutations = 10000) 

#Mantel statistic r: 0.4424 
#     Significance: 0.0009999 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.126 0.170 0.212 0.267 
#Permutation: free
#Number of permutations: 10000

#Hacer el Plot del IBD basado en Distancia Euclideana (Allsamples 16 islotes)
#https://stackoverflow.com/questions/54737007/how-do-i-fix-the-abline-warning-only-using-first-two-coefficients
#https://bookdown.org/hhwagner1/LandGenCourse_book/WE_5.html

plot(euclidean_distances_matrix,genetic_matrix, ylab="Fst",xlab="GeographicDistance(Euc)")
Dgeo=as.numeric(as.character(euclidean_distances_matrix))
Dgen=as.vector(genetic_matrix)
dens <- MASS::kde2d(Dgeo, Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
plot(Dgeo, Dgen, pch=20, cex=0.5,  
     xlab="Geographic Distance (Euc)", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(Dgen ~ as.numeric(as.character(Dgeo))))
lines(loess.smooth(Dgeo, Dgen), col="red")
title("IBD_Podarcis_lilfordi_Euc")

#Se realiza el TEST DE MANTEL usando ahora la distancia geodesica (Allsamples 16 islotes)

mantel.p_lilfordi_geo<-mantel(geo_distances_matrix,genetic_matrix,method="pearson",permutations=10000)
mantel.p_lilfordi_geo

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#mantel(xdis = geo_distances_matrix, ydis = genetic_matrix, method = "pearson",permutations = 10000) 

#Mantel statistic r: 0.4011 
 #     Significance: 0.0016998 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.129 0.175 0.219 0.274 
#Permutation: free
#Number of permutations: 10000

#Hacer el Plot del IBD basado en Distancia Geodesica (Allsamples 16 islotes)
plot(geo_distances_matrix,genetic_matrix, ylab="Fst",xlab="GeographicDistance(Geo)")
Dgeo=as.numeric(as.character(geo_distances_matrix))
Dgen=as.vector(genetic_matrix)
dens <- MASS::kde2d(Dgeo, Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
plot(Dgeo, Dgen, pch=20, cex=0.5,  
     xlab="Geographic Distance (Geo)", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(Dgen ~ as.numeric(as.character(Dgeo))))
lines(loess.smooth(Dgeo, Dgen), col="red")
title("IBD_Podarcis_lilfordi_Geo")

##### Vamos a correrlo entre Mallorca-Cabrera y Menorca aparte ########

############################## MENORCA TEST MANTEL (7 ISLOTES) #########################

library(geodist)
setwd("/Users/katherin_otalora/Documents/Scripts_R/FINAL_SCRIPT_ALLSAMPLES/MANTEL")
#coord_dec1 <- read.table("Coordenadas_decimales_Plilfordi_30_05_23.txt", header = TRUE)
coord_dec_menorca <- read.table("1Number_Coordenadas_decimales_Plilfordi_30_05_23_Menorca.txt", header = TRUE)
coord_dec_menorca

#coord_dec1
#euclidean distances in degrees:
euclidean_distances_menorca <- dist(coord_dec_menorca,diag=T, upper=T)
euclidean_distances_menorca
euclidean_distances_menorca_matrix  <-as.matrix(euclidean_distances_menorca, rownames=TRUE, rownames.value=NULL)
euclidean_distances_menorca_matrix 

#dist(coord_dec1, diag=T, upper=T)
#geodesic distances: 
geo_distances_matrix_menorca <-geodist(coord_dec_menorca)
geo_distances_matrix_menorca
#SALE UN WARNING PREGUNTAR??

setwd("/Users/katherin_otalora/Documents/Scripts_R/FINAL_SCRIPT_ALLSAMPLES/MANTEL")
#df1 <- read.table("FST_Plilfordi_30_05_23.txt", header = TRUE)
df_menorca <- read.table("1FST_Plilfordi_30_05_23_Menorca.txt", header = TRUE)
df_menorca
#df1
genetic_matrix_menorca <-as.matrix(df_menorca, rownames=TRUE, rownames.value=NULL)
#genetic_matrix <-as.matrix(df1, rownames=TRUE, rownames.value=NULL)

genetic_matrix_menorca

#Se realiza el TEST DE MANTEL usando primero la distancia euclideana (Menorca 7 islotes)

#install.packages ("vegan")
library (vegan)

mantel.p_lilfordi_euc_menorca<-mantel(euclidean_distances_menorca_matrix,genetic_matrix_menorca,method="pearson",permutations=10000)
mantel.p_lilfordi_euc_menorca

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#mantel(xdis = euclidean_distances_menorca_matrix, ydis = genetic_matrix_menorca,      method = "pearson", permutations = 10000) 

#Mantel statistic r: -0.4206 
#      Significance: 0.91468 

#Upper quantiles of permutations (null model):
 # 90%   95% 97.5%   99% 
#0.807 0.842 0.866 0.884 
#Permutation: free
#Number of permutations: 5039

#Hacer el Plot del IBD basado en Distancia Euclideana (Menorca 7 islotes)
#https://stackoverflow.com/questions/54737007/how-do-i-fix-the-abline-warning-only-using-first-two-coefficients
#https://bookdown.org/hhwagner1/LandGenCourse_book/WE_5.html

plot(euclidean_distances_menorca_matrix,genetic_matrix_menorca, ylab="Fst",xlab="GeographicDistance(Euc)")
Dgeo=as.numeric(as.character(euclidean_distances_menorca_matrix))
Dgen=as.vector(genetic_matrix_menorca)
dens <- MASS::kde2d(Dgeo, Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
plot(Dgeo, Dgen, pch=20, cex=0.5,  
     xlab="Geographic Distance (Euc)", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(Dgen ~ as.numeric(as.character(Dgeo))))
lines(loess.smooth(Dgeo, Dgen), col="red")
title("IBD_Podarcis_lilfordi_Menorca_Euc")

#Se realiza el TEST DE MANTEL usando ahora la distancia geodesica (Menorca 7 islotes)

mantel.p_lilfordi_geo_menorca<-mantel(geo_distances_matrix_menorca,genetic_matrix_menorca,method="pearson",permutations=10000)
mantel.p_lilfordi_geo_menorca

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#mantel(xdis = geo_distances_matrix_menorca, ydis = genetic_matrix_menorca,      method = "pearson", permutations = 10000) 

#Mantel statistic r: -0.4222 
 #     Significance: 0.91587 

 # 90%   95% 97.5%   99% 
#Upper quantiles of permutations (null model):
#0.808 0.843 0.868 0.885 
#Permutation: free
#Number of permutations: 5039

#Hacer el Plot del IBD basado en Distancia Geodesica (Allsamples 16 islotes)
plot(geo_distances_matrix_menorca,genetic_matrix_menorca, ylab="Fst",xlab="GeographicDistance(Geo)")
Dgeo=as.numeric(as.character(geo_distances_matrix_menorca))
Dgen=as.vector(genetic_matrix)
dens <- MASS::kde2d(Dgeo, Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
plot(Dgeo, Dgen, pch=20, cex=0.5,  
     xlab="Geographic Distance (Geo)", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(Dgen ~ as.numeric(as.character(Dgeo))))
lines(loess.smooth(Dgeo, Dgen), col="red")
title("IBD_Podarcis_lilfordi_Menorca_Geo")

############################## MALLORCA_CABRERA TEST MANTEL (9 ISLOTES) #########################

library(geodist)
setwd("/Users/katherin_otalora/Documents/Scripts_R/FINAL_SCRIPT_ALLSAMPLES/MANTEL")
#coord_dec1 <- read.table("Coordenadas_decimales_Plilfordi_30_05_23.txt", header = TRUE)
coord_dec_mallo_cab <- read.table("1Number_Coordenadas_decimales_Plilfordi_30_05_23_Mallorca_Cab.txt", header = TRUE)
coord_dec_mallo_cab

#coord_dec1
#euclidean distances in degrees:
euclidean_distances_mallo_cab <- dist(coord_dec_mallo_cab,diag=T, upper=T)
euclidean_distances_mallo_cab
euclidean_distances_mallo_cab_matrix  <-as.matrix(euclidean_distances_mallo_cab, rownames=TRUE, rownames.value=NULL)
euclidean_distances_mallo_cab_matrix 

#dist(coord_dec1, diag=T, upper=T)
#geodesic distances: 
geo_distances_matrix_mallo_cab <-geodist(coord_dec_mallo_cab)
geo_distances_matrix_mallo_cab
#SALE UN WARNING PREGUNTAR??

setwd("/Users/katherin_otalora/Documents/Scripts_R/FINAL_SCRIPT_ALLSAMPLES/MANTEL")
#df1 <- read.table("FST_Plilfordi_30_05_23.txt", header = TRUE)
df_mallo_cab <- read.table("1FST_Plilfordi_30_05_23_Mallorca_Cab.txt", header = TRUE)
df_mallo_cab
#df1
genetic_matrix_mallo_cab <-as.matrix(df_mallo_cab, rownames=TRUE, rownames.value=NULL)
#genetic_matrix <-as.matrix(df1, rownames=TRUE, rownames.value=NULL)

genetic_matrix_mallo_cab

#Se realiza el TEST DE MANTEL usando primero la distancia euclideana (Menorca 7 islotes)

#install.packages ("vegan")
library (vegan)

mantel.p_lilfordi_euc_mallo_cab<-mantel(euclidean_distances_mallo_cab_matrix,genetic_matrix_mallo_cab,method="pearson",permutations=10000)
mantel.p_lilfordi_euc_mallo_cab

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#mantel(xdis = euclidean_distances_mallo_cab_matrix, ydis = genetic_matrix_mallo_cab,      method = "pearson", permutations = 10000) 

#Mantel statistic r: -0.08175 
#      Significance: 0.59784 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.429 0.568 0.641 0.676 
#Permutation: free
#Number of permutations: 10000


#Hacer el Plot del IBD basado en Distancia Euclideana (Mallorca-Cab 9 islotes)
#https://stackoverflow.com/questions/54737007/how-do-i-fix-the-abline-warning-only-using-first-two-coefficients
#https://bookdown.org/hhwagner1/LandGenCourse_book/WE_5.html

plot(euclidean_distances_mallo_cab_matrix,genetic_matrix_mallo_cab, ylab="Fst",xlab="GeographicDistance(Euc)")
Dgeo=as.numeric(as.character(euclidean_distances_mallo_cab_matrix))
Dgen=as.vector(genetic_matrix_mallo_cab)
dens <- MASS::kde2d(Dgeo, Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
plot(Dgeo, Dgen, pch=20, cex=0.5,  
     xlab="Geographic Distance (Euc)", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(Dgen ~ as.numeric(as.character(Dgeo))))
lines(loess.smooth(Dgeo, Dgen), col="red")
title("IBD_Podarcis_lilfordi_Mallorca_Cabrera_Euc")

#Se realiza el TEST DE MANTEL usando ahora la distancia geodesica (Mallorca-Cabrera 9 islotes)

mantel.p_lilfordi_geo_mallo_cab<-mantel(geo_distances_matrix_mallo_cab,genetic_matrix_mallo_cab,method="pearson",permutations=10000)
mantel.p_lilfordi_geo_mallo_cab

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#mantel(xdis = geo_distances_matrix_mallo_cab, ydis = genetic_matrix_mallo_cab,      method = "pearson", permutations = 10000) 

#Mantel statistic r: -0.07399 
    #  Significance: 0.60034 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.409 0.552 0.625 0.658 
#Permutation: free
#Number of permutations: 10000


#Hacer el Plot del IBD basado en Distancia Geodesica (Allsamples 16 islotes)
plot(geo_distances_matrix_mallo_cab,genetic_matrix_mallo_cab, ylab="Fst",xlab="GeographicDistance(Geo)")
Dgeo=as.numeric(as.character(geo_distances_matrix_mallo_cab))
Dgen=as.vector(genetic_matrix_mallo_cab)
dens <- MASS::kde2d(Dgeo, Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
plot(Dgeo, Dgen, pch=20, cex=0.5,  
     xlab="Geographic Distance (Geo)", ylab="Genetic Distance")
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(Dgen ~ as.numeric(as.character(Dgeo))))
lines(loess.smooth(Dgeo, Dgen), col="red")
title("IBD_Podarcis_lilfordi_Mallorca_Cabrera_Geo")


#Fin
