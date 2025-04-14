######### ADMIXTURE ###########

########################################################################### DATA SET FROM GATK ###########################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################


### HACER UN .bed para ADMIXTURE en PLINK ######

cd /Users/katherin_otalora/Documents/Plink_1/plink

VCF=/Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf

#1. Descargas plink y en la misma carpeta pones tu vcf y corres la siguiente linea:

./plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 5 0.20 --out intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_prune_50_5_02

#2. PODAR EL ARCHIVO DE LOS SNPS CON LD:

./plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --extract intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_prune_50_5_02.prune.in --make-bed --pca --out /Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02

#3.TRANSFORMAR EL PRUNEDDATA ( ARCHIVO PODADO) A VCF DE NUEVO

./plink2  --bfile /Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02 --double-id --allow-extra-chr --set-missing-var-ids @:# --cow --recode vcf --out /Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02

grep -v "^#" /Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.vcf | wc -l

#Esto ya se habia hecho en el FILTRADO ver archivo de filtrado

################################################################### CORREMOS CON EL OUTPUT .bed ################################################################

#Correr Admixture ###### Es importante el .bim en la primera columna tenga ceros porque en ocasiones no corre ( esta columna al poner 0 dice que se desconoce el cromosoma)
#LOS DOCUMENTOS .BED .FAM .BIM deben estar en la misma carpeta que el ejecutable de admixture si no, no correra, si no hay que poner toda la linea de la ubicacion del ejecutable

xattr -cr ./admixture

cd /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0

./admixture

# Make the first column a bunch of zeros.

awk '{$1=0;print $0}' intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bim > intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bim.tmp
awk '{$1=0;print $0}' extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bim > extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bim.tmp
awk '{$1=0;print $0}' allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.bim > allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.bim.tmp
awk '{$1=0;print $0}' allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_LD_prunned_50_5_02.bim > allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_LD_prunned_50_5_02.bim.tmp

# Rename the .bim file

mv intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bim.tmp intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bim
mv extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bim.tmp extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bim
mv allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.bim.tmp allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.bim
mv allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_LD_prunned_50_5_02.bim.tmp allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_LD_prunned_50_5_02.bim

#Prueba para ver si corre bien admixture
./admixture intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bed 3

## run admixture for K 1-8 el Crossvalidation viene a 10 folds por defecto, nosotros hicimos pruebas con 50, 20, 10 y nos quedamos con 10
#CV=10 es cvfolds 10

#Para intSamples corremos con 8 poblaciones y cv=10
for K in 1 2 3 4 5 6 7 8; \
do /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/admixture --cv=10 /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/GATK/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bed   $K | tee /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/GATK/log${K}.out;
done

#Para extSamples corremos con 10 poblaciones y cv=10
for K in 1 2 3 4 5 6 7 8 9 10; \
do /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/admixture --cv=10 /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/GATK/extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.bed   $K | tee /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/GATK/log${K}.out;
done

#Para allSamples 0.90 corremos con 16 poblaciones y cv=10

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/admixture --cv=10 /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/GATK/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.bed   $K | tee /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/GATK/log${K}.out;
done

#Para allSamples 0.95 corremos con 16 poblaciones y cv=10

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/admixture --cv=10 /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/GATK/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_LD_prunned_50_5_02.bed   $K | tee /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/GATK/log${K}.out;
done

##################### R software intSample ###################################

#Se modifico el file 3.Q donde los primeros 4 individuos se colocaron cerca de cada una de las poblaciones a las que pertenecian, y ADRG se puso antes de AIRE

setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/intSamples_admixture")

tbl=read.table("intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.3_1.Q")
setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/intSamples_admixture")
name_plilfordi <- read.table(file = "popfile_1.txt", header = TRUE)

barplot(t(as.matrix(tbl)), col=rainbow(4),
               xlab="Individual #", ylab="Ancestry", border=NA, main="Admixture P.lilfordi")

               names.arg = "name_plilfordi")

#Se modifico el file 4.Q donde los primeros 4 individuos se colocaron cerca de cada una de las poblaciones a las que pertenecian, y ADRG se puso antes de AIRE.

tbl=read.table("intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.4_1.Q")
setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/intSamples_admixture")
name_plilfordi <- read.table(file = "popfile_1.txt", header = TRUE)

barplot(t(as.matrix(tbl)), col=rainbow(4),
                              xlab="Individual #", ylab="Ancestry", border=NA, main="Admixture P.lilfordi")

##################### R software extSample ###################################

#Se modifico el file 4.Q donde los primeros 4 individuos se colocaron cerca de cada una de las poblaciones a las que pertenecian, y ADRG se puso antes de AIRE

setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/extSamples_admixture")

tbl=read.table("extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.4_1.Q")
setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/extSample_admixture")
name_plilfordi <- read.table(file = "popfile_2.txt", header = TRUE)

barplot(t(as.matrix(tbl)), col=rainbow(4),
xlab="Individual #", ylab="Ancestry", border=NA, main="Admixture P.lilfordi")

#Se modifico el file 4.Q donde los primeros 4 individuos se colocaron cerca de cada una de las poblaciones a las que pertenecian, y ADRG se puso antes de AIRE

setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/extSamples_admixture")

tbl=read.table("extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.3_1.Q")
setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/extSamples_admixture")
name_plilfordi <- read.table(file = "popfile_2.txt", header = TRUE)

barplot(t(as.matrix(tbl)), col=rainbow(4),
xlab="Individual #", ylab="Ancestry", border=NA, main="Admixture P.lilfordi")

##################### R software allSample ###################################

#Se modifico el file 9.Q donde los primeros 4 individuos se colocaron cerca de cada una de las poblaciones a las que pertenecian, y ADRG se puso antes de AIRE

setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/allSamples_admixture_0.90")

tbl=read.table("allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.9_1.Q")
slist1 <- alignK(tbl[c(3,4,10)])
setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/allSamples_admixture_0.90")
name_plilfordi <- read.table(file = "popfile_3_nom.txt", header = TRUE)

barplot(t(as.matrix(tbl)), col=rainbow(9),
xlab="Individual #", ylab="Ancestry", border="grey50", main="Admixture P.lilfordi")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
# Load required libraries
library(ggplot2)

# Establecer el directorio de trabajo
setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/allSamples_admixture_0.90")

# Leer el archivo de resultados
tbl <- read.table("allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.9_1.Q")

# Leer el archivo de nombres de poblaciones
name_plilfordi <- read.table(file = "popfile_3_nom.txt", header = TRUE)

# Unir las columnas "CODE" y "POP" al dataframe tbl
tbl <- cbind(tbl, name_plilfordi)

# Create a long-format data frame for ggplot
df <- reshape2::melt(tbl, id.vars = c("CODE", "POP"))

# Create the bar plot
ggplot(df, aes(x = as.factor(variable), y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual #", y = "Ancestry", title = "Admixture P.lilfordi") +
  scale_fill_manual(values = rainbow(ncol(tbl))) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_grid(~CODE + POP, scales = "free_y", space = "free_y") +
  coord_flip()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

for (i in 1:length(name_plilfordi)){
  axis(1, at=median(which(sample_sites==name_plilfordi[i])), labels=name_plilfordi[i])}

barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

par(mar=c(10,4,4,4))

barplot(t(as.matrix(ordered[,2:9])), col=rainbow(6), border=NA,
        names.arg=barNaming(ordered$name_plilfordi), las=2)


#Se modifico el file 9.Q donde los primeros 4 individuos se colocaron cerca de cada una de las poblaciones a las que pertenecian, y ADRG se puso antes de AIRE
install.packages("tidyverse") 
library(tidyverse)
setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/allSamples_admixture_0.90")
tbl=read.table("allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.9.Q")

tbl = read.table("~/path-to/my_admixture_analysis.4.Q")
ord = tbl[order(tbl$V1,tbl$V2,tbl$V3,tbl$V4,tbl$V5,tbl$V6,tbl$V7,tbl$V8,tbl$V9),]
bp = barplot(t(as.matrix(ord)), 
             space = c(0.1),
             col=rainbow(9),
             xlab="Individual #", 
             ylab="Ancestry",
             border=NA)
samplelist <- read.csv("popfile_3.tx",
                       col_names = "sample")


###################################### OPCIONES NUEVAS PARA GRAFICAR POR ORDEN DE ASIGNACION #########################################################

library(tidyverse)

plot_data <- tbl %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V9) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))


ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  theme_classic()

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  facet_grid(~likely_assignment, scales = 'free', space = 'free')

######################################### 09 oct 2023 Graficar varios K con Admixture ##################

cd /Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/allSamples_admixture_0.90

#Output file
FILE=allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02

############### Prepare to plot in R #########################

#run it. in the terminal This will generate a tiff file ###
#Use this script. t requires four arguments, the prefix for the ADMIXTURE output files (-p ), the file with the species information (-i ), the maximum number of K to be plotted (-k 9), and a list with the populations or species separated by commas (-l <pop1,pop2...>). The list of populations provided with -l gives the order in which the populations or species shall be plotted. 

wget https://github.com/speciationgenomics/scripts/raw/master/plotADMIXTURE.r
chmod +x plotADMIXTURE.r

# the script plotADMIXTURE.r can be modified for better graphics

FILE=allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02
Rscript plotADMIXTURE.r -p $FILE -i popfile_3_1.txt -k 2 -l PORROS,AIRE,RAVL,PRCV,ADGR,REI,COLOM,DRAGONERA,CABRERA_H,CABRERA_L,COLOMER,ESCLTS,FORADADA,MOLTONA,GUARDIA,ESCURT

setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/allSamples_admixture_0.90")
plotADMIXTURE.r (qFile="allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.9.Q", metaData=popfile_3_1.txt, graph=NULL, palette="Dark2", colourN=8, structurePlot=TRUE, save=FALSE)

########################### POPHELPER #####################################

#found the answer, you need first to move all *Q files from each K into a single folder. Next, you should process this folder using 'pophelper' in R, using the following command:
  
remotes::install_github('royfrancis/pophelper')

library('pophelper')

setwd("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/allSamples_admixture_0.90/Q_data")

sfiles <- list.files(path=system.file("/Users/katherin_otalora/Documents/Admixture/dist/admixture_macosx-1.3.0/allSamples_admixture_0.90",package="pophelper"), full.names=T)
slist <- readQ (sfiles)

#Pophelper will create 'paramfiles' and combined txt files in folders per each K. Next, make sure you place clumpp.exe file into each folder, click on the executable file and it will generate three other files. Use the *-combined-merged, copy values and plot it in excel.

