###############################Correlacion entre variables Genetics y Ecologicas 2023 ################
#Script basado en: https://rpubs.com/osoramirez/316691

#Instalacion de programas
install.packages("PerformanceAnalytics")
install.packages("corrplot")
install.packages("visreg")
install.packages("ggplot2")

#LLamado de librerias
library(PerformanceAnalytics)
library(corrplot)
library(visreg)
library(ggplot2)

#Entramos a la carpeta donde tenemos nuestra informacion organizada.
setwd("/Users/katherin_otalora/Documents/Scripts_R/FINAL_SCRIPT_ALLSAMPLES/Correlacion_Regresion")
variables <- read.table("Tabla_size_ar_allsamples_07_06_23.txt", header = TRUE)
variables
#Volvemos nuestras varibles como numericas
Island_Area <- as.numeric(variables$Island_Area_ha)
Lizard_Density <- as.numeric(variables$Lizard_Density_ind_ha)
Maximum_Altitude <- as.numeric (variables$Maximum_Altitude_msnm)
Biotic_capacity <- as.numeric(variables$Index_Biotic_capacity)
Ar <- as.numeric(variables$Ar)
Pi <- as.numeric(variables$Pi)
H <- as.numeric(variables$Ho)

#Evaluamos una correlación rapida entre nuestras variables para ver el comportamiento entre las variables
#Cambiar el numero despues del : segun el numero de variables que tengamos en este caso 8
pairs.panels(variables[,2:8], scale=T, method="pearson", hist.col="orange")
pairs.panels(variables[,2:8], scale=T, method="spearman", hist.col="cyan")

###################### SHAPIRO.TEST PARAMETRICOS O NO PARAMETRICOS ###############################
#Veamos que variables son normales y cuales no, para poder transformarlas.
#Veamos primero si nuestrso datos son parametricos (distribuccion normal) o no parametricos
# Luego de saber si son o no parametricos podemos elegir si hacer Coor de Pearson (Parametricos) o Spearman (No parametricos)

shapiro.test(Island_Area) #No es normal p-value = 1.513e-06, hay que transformarla
shapiro.test(Lizard_Density) # Patron de normalidad p-value = 0.3685
shapiro.test(Maximum_Altitude) # No es normal p-value  p-value = 0.0001102, hay que transformarla
shapiro.test(Biotic_capacity) # Patron de normalidad  p-value = 0.3752
shapiro.test(Ar) # Patron de normalidad  p-value =0.2297
shapiro.test(Pi) # Patron de normalidad  p-value =0.3028
shapiro.test(H) # Patron de normalidad p-value = 0.1829

# Veamos si al transformar las variables ya siguen el patron de normalidad
shapiro.test(log(Island_Area)) # Patron de normalidad p-value = 0.5137
shapiro.test(log(Maximum_Altitude)) # Patron de normalidad p-value = 0.7595
shapiro.test(log(Lizard_Density)) # Patron de normalidad

################# EVALUACION DE CORRELACION Y REGRESION POR PARES ###############
######################### DATOS TRANSFORMADOS #############################
#Transformacion de mis datos con Log hacemos las correlacion para Island_Area y Ar
#Revisamos mediante un PLOT de correlacion si hay una relacion entre las variables
#Revisamos si hay una relacion entre las variables a la hora de transformarlas con log

################## log(Island_Area) #########################
# Nuestros datos es mejor tenerlos en un data.frame

######## Island area and Ar ############
dat1 <- data.frame(log(Island_Area),Ar)
chart.Correlation(dat1)
cor.test(log(Island_Area),Ar)  # cor.test analiza la significancia de la correlacion, utilizando el contraste
#0.7999091, p-value=0.0001992 
M1 <- cor(dat1)
corrplot(M1, method = "ellipse") #Para graficarlo para un paper
0.7999091^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#  0.6398546 ##En nuestros datos el coeficiente de determinación es  0.6398546 = 0,63 o 63%, de lo que se deduce que aproximadamente el 37% del sizeisland no esta relacionado con el Ar (diversidad genetica).

######## Island area and Pi ############
dat2 <- data.frame(log(Island_Area),Pi)
chart.Correlation(dat2)
cor.test(log(Island_Area),Pi)  # cor.test analiza la significancia de la correlacion, utilizando el contraste
#0.3578676, p-value = 0.1735
M2 <- cor(dat2)
corrplot(M2, method = "ellipse") #Para graficarlo para un paper
0.3578676^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#  0.1280692 ##En nuestros datos el coeficiente de determinación es  0.1280692 = 0,12 o 12%, de lo que se deduce que aproximadamente el 88% del sizeisland no esta relacionado con el Ar (diversidad genetica).

######## Island area and H ############
dat3 <- data.frame(log(Island_Area),H)
chart.Correlation(dat3)
cor.test(log(Island_Area),H)  # cor.test analiza la significancia de la correlacion, utilizando el contraste
#0.3435916, p-value = 0.1926
M3 <- cor(dat3)
corrplot(M3, method = "ellipse") #Para graficarlo para un paper
0.3435916^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#  0.1180552 ##En nuestros datos el coeficiente de determinación es  0.1180552 = 0,11 o 11%, de lo que se deduce que aproximadamente el 89% del sizeisland no esta relacionado con el Ar (diversidad genetica).

################## log(Lizard_Density) #########################

######## Lizard_Density and Ar ############ Usamos SPEARMAN (Datos no normales)
dat4 <- data.frame(log(Lizard_Density),Ar)
chart.Correlation(dat4)
cor.test(log(Lizard_Density),Ar) #Se hizo Pearson para ver si hay diferencias con Spearman
#0.1038799, p-value = 0.7018 r2=0.01079103
cor.test(log(Lizard_Density),Ar, method = "spearman")
#-0.1030929, p-value = 0.704 r2= -0.01062815
M4 <- cor(dat4)
corrplot(M4, method = "ellipse") #Para graficarlo para un paper
0.1038799^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#  0.01079103 ##En nuestros datos el coeficiente de determinación es  0.05552663 = 0,05 o 5%, de lo que se deduce que aproximadamente el 95% del density no esta relacionado con el Pi (diversidad genetica).

######## Lizard_Density and Pi ############ 
dat5 <- data.frame(log(Lizard_Density),Pi)
chart.Correlation(dat5)
cor.test(log(Lizard_Density),Pi) #Se hizo Pearson para ver si hay diferencias con Spearman
#0.002997628 , p-value = 0.9912 r2=8.985774e-06
cor.test(log(Lizard_Density),Pi, method = "spearman")  # cor.test analiza la significancia de la correlacion, utilizando el contraste
#-0.07069227, p-value = 0.7947 r2=-0.004997397
M5 <- cor(dat5)
corrplot(M5, method = "ellipse") #Para graficarlo para un paper
0.002997628^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#  8.985774e-06 #En nuestros datos el coeficiente de determinación es  0.08503383 = 0,08 o 8%, de lo que se deduce que aproximadamente el 92% del density no esta relacionado con el Pi(diversidad genetica).

######## Lizard_Density and H ############ 
dat6 <- data.frame(log(Lizard_Density),H)
chart.Correlation(dat6)
cor.test(log(Lizard_Density),H) #Se hizo Pearson para ver si hay diferencias con Spearman
#0.02696819  , p-value = 0.921, r2=0.0007272833
cor.test(log(Lizard_Density),H, method = "spearman")  # cor.test analiza la significancia de la correlacion, utilizando el contraste
#-0.2106041,p-value = 0.4337, r2=-0.04435409
M6 <- cor(dat6)
corrplot(M6, method = "ellipse") #Para graficarlo para un paper
0.02696819^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#  0.0007272833 ##En nuestros datos el coeficiente de determinación es  0.09115261 = 0,09 o 9%, de lo que se deduce que aproximadamente el 91% del density no esta relacionado con el Pi(diversidad genetica).

############################ log(Maximum_Altitude)#######################

######## Maximum_Altitude and Ar ############ 
dat7 <- data.frame(log(Maximum_Altitude),Ar)
chart.Correlation(dat7)
cor.test(log(Maximum_Altitude),Ar)  # cor.test analiza la significancia de la correlacion d Pearson, utilizando el contraste
#0.6276103 , p-value = 0.009247
M7 <- cor(dat7)
corrplot(M7, method = "ellipse") #Para graficarlo para un paper
0.6276103^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#0.3938947 #En nuestros datos el coeficiente de determinación es  0.3938947 = 0,39 o 39%, de lo que se deduce que aproximadamente el 61% del sizeisland esta no relacionado con el Ar (diversidad genetica).

######## Maximum_Altitude and Pi ############ 
dat8 <- data.frame(log(Maximum_Altitude),Pi)
chart.Correlation(dat8)
cor.test(log(Maximum_Altitude),Pi)  # cor.test analiza la significancia de la correlacion d Pearson, utilizando el contraste
#-0.0273674  , p-value = 0.9199
M8 <- cor(dat8)
corrplot(M8, method = "ellipse") #Para graficarlo para un paper
-0.0273674^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#-0.0007489746 #En nuestros datos el coeficiente de determinación es  -0.0007489746 = 0% , de lo que se deduce que aproximadamente el 99%%  no relacionado con el Pi (diversidad genetica).

######## Maximum_Altitude and H ############ 
dat9 <- data.frame(log(Maximum_Altitude),H)
chart.Correlation(dat9)
cor.test(log(Maximum_Altitude),H)  # cor.test analiza la significancia de la correlacion d Pearson, utilizando el contraste
#0.02748883  , p-value = 0.9195
M9 <- cor(dat9)
corrplot(M9, method = "ellipse") #Para graficarlo para un paper
0.02748883^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#0.0007556358 #En nuestros datos el coeficiente de determinación es 0.0007556358 = 0,1% , de lo que se deduce que aproximadamente el 99%%  no relacionado con el H (diversidad genetica).

############################  Biotic_capacity ############################# 

########## Biotic_capacity and Ar############################# 
dat10 <- data.frame(Biotic_capacity,Ar)
chart.Correlation(dat10)
cor.test(Biotic_capacity,Ar)  # cor.test analiza la significancia de la correlacion d Pearson, utilizando el contraste
#0.7605487  , p-value = 0.0006251
M10 <- cor(dat10)
corrplot(M10, method = "ellipse") #Para graficarlo para un paper
0.7605487^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#0.5784343 #En nuestros datos el coeficiente de determinación es  0.5784343 = 0,57% , de lo que se deduce que aproximadamente el 57%% esta relacionado con el Pi (diversidad genetica).

########## Biotic_capacity and Pi############################# 
dat11 <- data.frame(Biotic_capacity,Pi)
chart.Correlation(dat11)
cor.test(Biotic_capacity,Pi)  # cor.test analiza la significancia de la correlacion d Pearson, utilizando el contraste
#0.2188229  , p-value = 0.4155
M11 <- cor(dat11)
corrplot(M11, method = "ellipse") #Para graficarlo para un paper
0.2188229^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#0.04788346 #En nuestros datos el coeficiente de determinación es  0.04788346 = 4% , de lo que se deduce que aproximadamente el 96%% NO esta relacionado con el Pi (diversidad genetica).

########## Biotic_capacity and H############################# 
dat12 <- data.frame(Biotic_capacity,H)
chart.Correlation(dat12)
cor.test(Biotic_capacity,H)  # cor.test analiza la significancia de la correlacion d Pearson, utilizando el contraste
#0.2852421  , p-value = 0.2842
M12 <- cor(dat12)
corrplot(M12, method = "ellipse") #Para graficarlo para un paper
0.2852421^2 #Calculamos EL COEFICIENTE DE DETERMINACION R2
#0.08136306 #En nuestros datos el coeficiente de determinación es  0.08136306 = 8% , de lo que se deduce que aproximadamente el 92%% NO esta relacionado con el H (diversidad genetica).


#Otra forma de graficar un plot de correlacion entre las variables transformadas
#pairs(log(Island_Area)  ~ Ar)

############################## REGRESION LINEAL SIMPLE Log(Island_Area) y Ar ########################################
#La regresión supone que hay una variable fija, controlada por el investigador (es la variable independiente o predictora), y otra que no está controlada (variable respuesta o dependiente). 
#La correlación supone que ninguna es fija: las dos variables están fuera del control de investigador.

########Se hace la regresion lineal simple siempre con los datos de Island area transformados ###

ej1 <- data.frame (Ar, (log (Island_Area)))
reg1 <- lm(Ar ~ (log (Island_Area)), data = ej1)
summary(reg1)

#La ecuacion de prediccion obtenida por la regresion seria:
Ar = 1.341878 + 0.014147  * Island_Area

#Se aprueba con un p-value: 0.0001992
anova(reg1)

#Observar si el anova es significativo. 
#Muy importante tener este resultado en las pruebas de la regresión, dado que esto se ajusta a lo que buscamos en la hipótesis de la regresión.
#Ho = pendientes es igual a cero, H1 = pendientes es diferente a cero

#Verificacion de Supuestos

par(mfrow = c(2, 2))
plot(reg1)

#Prueba de normalidad de los residuos 
shapiro.test(resid(reg1)) #Residuos normales  p-value = 0.8073

#Forma grafica con ggplot2
library(ggplot2)#ggplot2 es una extension poderosa para graficar
ggplot(ej1, aes(x=(log (Island_Area)), y=Ar)) +
  geom_point(shape=1) +    # genera circulos en el grafico
  geom_smooth(method=lm)   # adjunta la linea de regresion por defecto es al 95% de confianza

############################## REGRESION LINEAL SIMPLE Log (Maximum_Altitude) and Ar ########################################

#La regresión supone que hay una variable fija, controlada por el investigador (es la variable independiente o predictora), y otra que no está controlada (variable respuesta o dependiente). La correlación supone que ninguna es fija: las dos variables están fuera del control de investigador.

########Se hace la regresion lineal simple siempre con los datos de Maximum_Altitude transformados ###

ej2 <- data.frame (Ar, (log (Maximum_Altitude)))
reg2 <- lm(Ar ~ (log (Maximum_Altitude)), data = ej2)
summary(reg2)

#La ecuacion de prediccion obtenida por la regresion seria:
Ar = 1.307589 + 0.020617  * Maximum_Altitude

#Se aprueba con un p-value: 0.009247
anova(reg2)

#Observar si el anova es significativo. 
#Muy importante tener este resultado en las pruebas de la regresión, dado que esto se ajusta a lo que buscamos en la hipótesis de la regresión.
#Ho = pendientes es igual a cero, H1 = pendientes es diferente a cero

#Verificacion de Supuestos

par(mfrow = c(2, 2))
plot(reg2)

#Prueba de normalidad de los residuos 
shapiro.test(resid(reg2)) #Residuos normales  p-value = 0.8073

#Forma grafica con ggplot2
library(ggplot2)#ggplot2 es una extension poderosa para graficar
ggplot(ej2, aes(x=(log (Maximum_Altitude)), y=Ar)) +
  geom_point(shape=1) +    # genera circulos en el grafico
  geom_smooth(method=lm)   # adjunta la linea de regresion por defecto es al 95% de confianza

############################## REGRESION LINEAL SIMPLE Biotic_capacity and Ar ########################################

#La regresión supone que hay una variable fija, controlada por el investigador (es la variable independiente o predictora), y otra que no está controlada (variable respuesta o dependiente). La correlación supone que ninguna es fija: las dos variables están fuera del control de investigador.

########Se hace la regresion lineal simple siempre con los datos de Maximum_Altitude transformados ###

ej3 <- data.frame (Ar, Biotic_capacity)
reg3 <- lm(Ar ~ Biotic_capacity, data = ej3)
summary(reg3)

#La ecuacion de prediccion obtenida por la regresion seria:
Ar = 1.329359 + 0.008697  * Biotic_capacity

#Se aprueba con un p-value: 0.009247
anova(reg3)

#Observar si el anova es significativo. 
#Muy importante tener este resultado en las pruebas de la regresión, dado que esto se ajusta a lo que buscamos en la hipótesis de la regresión.
#Ho = pendientes es igual a cero, H1 = pendientes es diferente a cero

#Verificacion de Supuestos

par(mfrow = c(2, 2))
plot(reg3)

#Prueba de normalidad de los residuos 
shapiro.test(resid(reg3)) #Residuos normales  p-value = 0.8073

#Forma grafica con ggplot2
library(ggplot2)#ggplot2 es una extension poderosa para graficar
ggplot(ej3, aes(x=(log (Biotic_capacity)), y=Ar)) +
  geom_point(shape=1) +    # genera circulos en el grafico
  geom_smooth(method=lm)   # adjunta la linea de regresion por defecto es al 95% de confianza


#################### Regresion lineal multiple (Backward) Ar, Island_Area ,Biotic_capacity and log(Maximum_Altitude)) NEW DATA ##################################################

# Generamos un modelo de regresion multiple (variables independientes) 
#Solo le ponemos log a Island Area porque a Lizard Density no hace falta.porque sino (y=a+^b)
ej4 <- data.frame(log(Island_Area), Ar, Biotic_capacity, log(Maximum_Altitude))
reg4 <- lm(Ar ~ log(Island_Area) + Biotic_capacity + log(Maximum_Altitude) , data = ej4)
summary(reg4)
ej5 <- data.frame( log(Island_Area), Ar, Biotic_capacity)
reg5 <- lm(Ar ~ log(Island_Area) + Biotic_capacity, data = ej5)
summary(reg5)
ej6 <- data.frame( log(Island_Area), Ar, log(Maximum_Altitude))
reg6 <- lm(Ar ~ log(Island_Area) + log(Maximum_Altitude), data = ej6)
summary(reg6)
ej7 <- data.frame( log(Island_Area), Ar)
reg7 <- lm(Ar ~ log(Island_Area), data = ej7)
summary(reg7)

# Revisamos el modelo
anova(reg4)
anova(reg5)
anova(reg6)
anova(reg7)


# Verificamos los supuestos
par(mfrow = c(2, 2))
plot(reg2)

# Revisamos la distribucion de los residuos
shapiro.test(residuals(reg2))
