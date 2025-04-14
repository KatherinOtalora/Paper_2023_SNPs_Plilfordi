######################################## TREEMIX ###############################################

cd /Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13
./configure
make
make install

###### Build the ML graph without root with migration (-m) 16 and bootstrap (-500 genera replicas de bootstrap remuestreando bloques de 500) uniendo CABRERA
#Tomamos el output que saca STACKs y de ahi le ponemos .txt y luego hacemos lo siguientes en la terminal:
# Se volvío a correr STACKs y ahi se unio cabrera para que diera el output de treemix, se le puso .txt y se comprimio en .gz para que lo lea treemix 
#Importante quitarle la primera linea de informacion del .txt sobre stacks porque el gzip no la reconoce bien y da error

gzip allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_1_cabrera.p.txt

#Vamos a la carpeta 
cd /Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/INPUT_Treemix_allSamples_0.90_LD_Script


#################################JULIO DE 2023##########################################

###### Build the ML graph without root with migration (-m) 5 and bootstrap (-500 genera replicas de bootstrap remuestreando bloques de 500)
#Corremos el siguiente loop sugerido por Ling Hu et al., y curso de speciation genomics http://speciationgenomics.github.io


cd /Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/25_JUL_23

for rep in {1..10}
 do
  for m in {1..5}
   do
        treemix -i allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_1_cabrera.p.txt.gz -m $m -bootstrap -k500 -global -se -o out.stem.allsamples_0.90_LD_withoutroot.m16.m10.global.se.B500.${rep}.${m}
    done
done


#Vamos a R y ponemos en la misma carpeta lo que queremos graficar con los resultados
#El paquete de plotting funcs con los datos de treemix en la misma carpeta y correr el siguiente script:

install.packages("R.utils")
library(RColorBrewer)
library(R.utils)
source("/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/src/plotting_funcs.R")

setwd ("/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/25_JUL_23")
pplot_tree("/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/src/out.stem.allsamples_0.90_LD_withoutroot.m16.m10.global.se.B500.10.1")
#Se selecciona el 10.1 porque según OptM la probabilidad de migración seria de un evento y en se debe seleccionar las ultimas iteraciones y con mayor numero de "m weight", que correspondería a este resultado con  0.3093
#El nombre dice m16 pero fueron de 1 a 5 migraciones con 10 iteraciones

...............................................F3 Y F4 statistics...........................................................................

threepop -i allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_1_cabrera.p.txt.gz -k 500 > threepop_allSamples_0.90_LD_JUL_23.txt

fourpop -i allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_1_cabrera.p.txt.gz -k 500 > fourpop_allSamples_0.90_LD_JUL_23.txt


############################ 28 AGOSTO 2023 -K500 con root Dragonera ######################################

###### Build the ML graph without root with migration (-m) 5 and bootstrap (-1000 genera replicas de bootstrap remuestreando bloques de 1000)


cd /Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/28_AGO_23

for rep in {1..10}
 do
  for m in {0..5}
   do
        treemix -i allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_1_cabrera.p.txt.gz -m $m -bootstrap -k500 -root DRAGONERA -global -se -o out.stem.allsamples_0.90_LD_withoutroot.m5.r10.global.se.B500.${rep}.${m}
    done
done

#Vamos a R y ponemos en la misma carpeta lo que queremos graficar con los resultados
#El paquete de plotting funcs con los datos de treemix en la misma carpeta y correr el siguiente script:

install.packages("R.utils")
library(RColorBrewer)
library(R.utils)
source("/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/src/plotting_funcs.R")

setwd ("/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/28_AGO_23")
plot_tree("/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/src/out.stem.allsamples_0.90_LD_withoutroot.m5.r10.global.se.B500.10.1")

#No se vuelve a correr el threepop y fourpop porque no cambia en nada.

