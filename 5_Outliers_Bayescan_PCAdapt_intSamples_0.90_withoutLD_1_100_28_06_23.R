##################### Outliers WITHOUTLD INTSAMPLES K=8 #############################

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

cd /Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/intSamples_0.90_withoutLD_bayescan_1_100_K=8

#Corremos el Bayescan con 1:100
#Testeamos: 20 short pilots runs with 5000 integration, burn in set to 5 x 104 and thinning interval 10 and prior odds of 100.

/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/BayeScan2.1_macos64bits /Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/intSamples_filtredwithoutLD_bayescan_1.txt -od ./ -threads 8 \ -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100

#La primera línea consiste en las opciones básicas para asegurar que el programa se ejecuta correctamente - es decir, dónde está el archivo de entrada, en qué directorio escribir la salida (./ significa el directorio actual) y el número de hilos o núcleos que podemos utilizar (4 en este caso).

#La segunda línea especifica las opciones para los análisis y la cadena de Markov Monte Carlo y son las siguientes

#-n es el número de iteraciones que queremos muestrear del MCMC (5000 en este caso)
#-thin es el intervalo de thin, fijado en 10; esto significa que dejamos que el MCMC se ejecute durante 10 pasos entre cada muestra.
#-nbp es el número de ejecuciones piloto que Bayescan ejecuta para ajustar los parámetros del MCMC.
#-pilot es la duración de cada una de las ejecuciones piloto, 5000 iteraciones aquí.
#-burn es la longitud del burnin (es decir, la parte del MCMC que se descarta), aquí es de 50 000
#prior odds es un valor importante porque tiene un efecto sobre los falsos positivos: usamos 100 para ser conservadores, 10 es menos conservador

#Luego de correr saldran varios outputs en la misma carpeta donde se ubico binaries.

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

setwd("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/intSamples_withoutLD_1_100_K=8")

bayescan=read.table("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/intSamples_withoutLD_1_100_K=8/intSamples_filtredwithoutLD_bayescan_1_fst.txt")

#La primera columna del marco de datos de bayescan es el ID del SNP.
#Las tres columnas siguientes (prob, log10(P0), y qval) están relacionadas con la prueba de adaptación local considerando el logaritmo de las probabilidades posteriores - log10(PO) - y el valor q para el modelo con selección.
#La quinta columna da el tamaño del efecto específico del locus (parámetro alfa). La última proporciona el FST específico del locus promediado en todas las poblaciones.

#Descargue la lista de SNPs en el orden correcto. El formato .geste tiene los SNPs en el mismo orden que el vcf utilizado para producir el formato .geste. Por lo tanto, puede utilizar este comando en bash para extraer la tercera columna que contiene la información de identificación de cada SNP en su vcf.
#Importante para el archivo orginal de vcf hay que cortar la columna 1 y 2 y luego juntarlas porque en la columna tres no hay informacion solo un . (Se juntan en algun editor de texto con _ entre medio de la culmna 1 y 2)

grep -v "#" intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf | cut -f 1,2 > intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_snps_1_100_K=8.txt

#A continuación, importe la información de los SNPs.

SNPb=read.table("intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_snps_1_100_K=8.txt", header=FALSE)

#Combinar los nombres de los valores atípicos con los resultados del marco de datos del bayescan.

bayescan=cbind(SNPb, bayescan)

#Cambiar el nombre de las columnas.

colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST")

#Escribe los resultados.

write.table(bayescan, "intSamples_0.90_filtred_withoutLD_bayescan_results_1_100_K=8.txt", quote=FALSE, sep="\t", row.names=FALSE)

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
#       70851         1236        71185

#Escriba los resultados de los SNPs potencialmente bajo selección (valor q < 0,01).

write.table(neutral, "intSamples_1_100_K=8_neutral.txt", row.names=F, quote=F)
write.table(balancing, "intSamples_1_100_K=8_balancing.txt", row.names=F, quote=F)
write.table(positive, "intSamples_1_100_K=8_positive.txt", row.names=F, quote=F)

#Transformación Log del valor Q para crear el gráfico ggplot.

range(bayescan$Q_VALUE)

#[1] 0.0001 0.9728

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

plot_bayescan("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/intSamples_withoutLD_1_100_K=8/intSamples_filtredwithoutLD_bayescan_1_fst.txt")

$outliers
   [1]      6     16     18     23     26     27     30     59     61     98    120    179    208    209    219
  [16]    282    283    379    396    398    666    667    668    994   1469   1512   1594   1771   1781   1814
  [31]   1840   1885   2381   2395   2416   2645   2655   3249   3276   3592   3711   3807   3810   3813   3842
  [46]   4007   4039   4361   4951   4953   5047   5052   5578   5599   5600   5601   5605   5610   5622   5624
  [61]   5748   5839   5940   5943   5944   5945   6124   6273   6463   6497   6501   6514   6537   6585   6639
  [76]   6799   6812   6976   7181   7409   7410   7411   7412   7555   7607   7609   7613   7639   7656   7905
  [91]   7912   7944   8407   8491   8505   8539   8540   8541   8969   8970   8987   8988   9007   9390   9391
 [106]   9419   9685   9837   9924  10023  10072  10119  10120  10189  10192  10318  10367  10475  10576  10708
 [121]  10711  10718  11066  11396  11552  11985  12033  12034  12219  12433  12434  12435  12437  12439  12440
 [136]  12441  12442  12452  12453  12465  12470  12480  12486  12487  12498  12504  12510  12512  12515  12517
 [151]  12518  12519  12535  12542  12543  12572  12598  12831  13096  13203  13288  13416  13448  13459  13460
 [166]  13469  13472  13475  13477  13478  13483  13484  13538  13546  13563  13564  13569  13586  13587  13600
 [181]  13634  13645  13701  13703  13708  13746  13747  13849  14013  14015  14016  14042  14054  14121  14232
 [196]  14233  14273  14281  14282  14612  14694  14711  14715  14716  14754  14784  14813  14952  15072  15073
 [211]  15135  15342  15411  15460  15606  15712  15776  16084  16180  16223  16292  16294  16662  16664  16728
 [226]  16742  16792  16911  17013  17114  17324  17412  17771  17879  18041  18253  18275  18277  18279  18359
 [241]  19084  19141  19142  19143  19158  19159  19261  19264  19265  19279  19379  19471  19708  19753  19755
 [256]  19760  19806  19866  19952  20472  20473  20791  21141  21142  21169  21290  21410  21749  21770  21953
 [271]  21986  22046  22197  22336  22355  22484  22493  22596  22602  22699  22929  22998  22999  23031  23367
 [286]  23404  23412  23423  23424  23431  23432  23433  23435  23436  23437  23438  23439  23440  23441  23442
 [301]  23443  23444  23445  23446  23447  23448  23449  23470  23471  23474  23491  23492  23499  23503  23506
 [316]  23521  23522  23523  23525  23527  23529  23530  23531  23533  23536  23537  23538  23540  23541  23542
 [331]  23545  23547  23548  23550  23551  23554  23557  23561  23568  23585  23588  23592  23593  23594  23596
 [346]  23607  23619  23636  23645  23665  23704  23708  23713  23774  23839  23849  23869  23871  23985  23993
 [361]  23994  24001  24108  24111  24113  24114  24199  24285  24297  24475  24476  24509  24512  24513  24775
 [376]  24827  24836  24855  24878  24907  25278  25944  25951  26568  26599  26725  26729  26730  26755  26825
 [391]  26873  26908  26910  26911  27406  27507  27681  27700  27788  28024  28025  28033  28088  28486  28680
 [406]  28813  28908  29027  29054  29876  29993  30049  30175  30489  30614  30645  30700  30702  30719  30765
 [421]  30766  30983  30986  31100  31187  31297  31743  31965  32165  32349  32379  32404  32550  32606  32627
 [436]  32629  32630  33064  33191  33268  33315  33365  33570  33664  33825  33861  33936  33941  34033  34647
 [451]  34805  34950  34952  35044  35168  35170  35184  35204  35209  35300  35326  35481  35930  36033  36059
 [466]  36166  36204  36284  36574  36705  36745  36784  36820  37111  37263  37348  37557  38028  38207  39083
 [481]  39345  39373  39667  39999  40060  40266  40324  40327  40430  40438  40488  40506  40569  40790  40793
 [496]  40823  40883  40897  40904  41019  41021  41226  41377  41532  41643  41669  41713  41955  42133  42152
 [511]  42260  42384  42476  42577  42764  42765  42808  43337  43538  43586  43587  43592  43600  43601  43614
 [526]  43615  43616  43617  43630  43661  43884  44277  44316  44480  44501  44606  44652  44667  44810  45045
 [541]  45050  45053  45054  45755  46163  46941  47090  47157  47476  47500  47501  47651  47754  48038  48523
 [556]  48717  48755  48802  48946  49239  49285  49371  49536  49637  50232  50339  50360  50375  50456  50467
 [571]  50828  50832  51046  51170  51217  51355  51391  51526  51575  51855  51875  51876  51909  51917  51918
 [586]  52136  52417  52564  52567  52572  52662  52835  53005  53073  53084  53086  53091  53101  53134  53140
 [601]  53149  53159  53185  53187  53193  53213  53253  53302  53305  53313  53448  53593  53875  53897  54023
 [616]  54043  54166  54466  54533  54725  54846  54869  54889  54970  55036  55117  55199  55273  55276  55430
 [631]  55839  55887  56460  56647  56949  56951  57304  57312  57569  57579  58109  58129  58134  58664  58735
 [646]  58736  58744  58939  59614  59643  59809  59859  59860  59861  60444  60513  60948  61494  61496  61534
 [661]  61545  61546  61548  61563  61564  61742  61771  62049  62463  62718  62731  62827  63434  63556  63557
 [676]  63563  64074  64278  64474  64828  65197  65200  65844  66163  66165  66201  66202  66252  66253  66445
 [691]  66471  66508  66668  66669  66670  66671  67123  67164  67323  67420  67452  67914  67928  68112  68136
 [706]  68252  68299  68793  69137  69339  69470  69508  69511  69790  70079  70149  70396  70491  70493  70500
 [721]  70733  70776  70812  70817  70821  70822  70823  70827  70830  70831  70834  70840  70856  70939  71373
 [736]  71386  71558  71559  71613  71616  71718  71719  71720  71762  71914  71973  72126  72152  72281  72699
 [751]  72909  72911  73096  73097  73328  73642  74065  74067  74069  74123  74278  74279  74298  74355  74356
 [766]  74377  74449  74468  74469  74806  74831  75547  75686  75956  76075  76172  76320  76402  76411  76461
 [781]  76642  76773  76780  76921  76927  76930  76950  77048  77207  77236  77564  77574  77579  77637  77664
 [796]  77667  77808  77930  77936  77971  78254  78283  78340  78444  78445  78707  79185  79220  79544  79563
 [811]  79564  79595  79774  80142  80181  80182  80771  80772  80777  80778  80780  81051  81336  81337  81343
 [826]  81349  81350  81359  81360  81365  81367  81368  81369  81370  81384  81398  81411  81412  81418  81439
 [841]  81494  81495  81871  81910  82051  82057  82062  82068  82078  82079  82080  82085  82097  82100  82102
 [856]  82103  82105  82115  82118  82225  82731  82740  82998  83065  83650  83807  84124  84132  84205  84251
 [871]  84491  84804  84967  84980  85097  85481  85512  85768  85836  86208  87166  87267  87592  88006  88092
 [886]  88159  88192  88229  88264  88299  88331  88339  88538  88615  88659  88660  88662  88663  88668  88670
 [901]  88674  88678  88683  88696  88709  88734  88740  88749  88790  88803  88831  88853  88886  89512  90062
 [916]  90063  90065  90066  90732  90946  91041  91347  91484  91487  91488  91497  91539  91603  91688  91689
 [931]  91707  91919  91977  92023  92314  92420  92550  92691  92735  92867  92869  92871  92895  93278  93308
 [946]  93320  93426  93433  93532  93534  93677  93725  93811  94321  94323  94324  94360  94648  94649  94723
 [961]  95030  95244  95326  95476  95653  95656  95657  95658  95661  95662  95663  95666  95672  95675  95682
 [976]  95683  95686  95688  95689  95690  95708  95709  95710  95723  95724  95726  95727  95730  95734  95757
 [991]  95761  95768  95913  96218  96505  96610  96644  97439  97616  97617  97621  97623  97624  97625  97801
[1006]  97848  97850  97977  98041  98053  98055  98073  98155  98157  98181  98514  98723  98837  98839  98842
[1021]  99126  99255  99671  99797  99822  99940 100092 100126 100142 100155 100189 100339 100537 100570 100707
[1036] 100771 100838 100911 101153 101155 101156 101157 101158 101160 101161 101163 101165 101168 101170 101185
[1051] 101199 101328 101404 101409 101640 101642 101658 101682 102009 102148 102393 102603 103046 103300 103726
[1066] 103796 103817 104201 104206 104508 104509 104773 104882 104984 105114 105149 105258 105306 105438 105628
[1081] 105711 106071 106239 106564 106893 106896 106897 106898 106900 106904 106908 106910 106911 106912 106944
[1096] 106945 106946 106956 106976 106977 107004 107054 107055 107057 107060 107126 107668 107761 107766 107844
[1111] 107874 107891 108107 108122 108343 108550 108614 108711 108713 108714 108715 108858 108890 108983 109088
[1126] 109307 109463 109467 109497 109498 109610 110606 110899 111021 111362 111402 111465 111466 111470 111471
[1141] 111472 111475 111488 111540 111541 111542 111543 111547 111548 111573 111574 111577 111578 111599 111602
[1156] 111633 111822 111823 111824 111825 111955 112004 112567 112947 113006 113075 113116 113117 113118 113120
[1171] 113126 113131 113137 113141 113142 113152 113174 113176 113183 113535 113583 113849 113910 114015 114443
[1186] 114755 114797 114802 114985 114999 115420 115735 115740 115743 115931 115935 115937 115998 116016 116507
[1201] 116688 116758 116768 116911 117199 117200 117321 117375 117438 117439 117514 117670 117671 117718 117852
[1216] 117864 118247 118499 118573 118574 118955 119035 119664 119947 120003 120051 120066 120072 120080 120085
[1231] 120086 120087 120091 120093 120127 120164 120199 120210 120222 120226 120236 120237 120238 120239 120240
[1246] 120241 120242 120243 120244 120245 120246 120248 120249 120250 120251 120252 120253 120254 120255 120256
[1261] 120257 120258 120260 120261 120269 120270 120272 120278 120279 120282 120288 120290 120291 120294 120296
[1276] 120317 120323 120325 120329 120335 120337 120338 120341 120342 120362 120374 120375 120377 120390 120460
[1291] 120462 120465 120469 120484 120507 120514 120524 120526 120532 120566 120613 120615 120637 120640 120641
[1306] 120668 120670 120672 120687 120688 120755 120756 120762 120763 120764 120772 120778 120782 120843 120845
[1321] 120869 120899 120982 121024 121121 121129 121130 121160 121227 121228 121344 121444 121445 121447 121464
[1336] 121465 121499 121500 121502 121508 121514 121517 121547 121556 121615 121640 121641 121711 121712 121730
[1351] 121735 121746 121750 121766 121785 121808 121813 121847 121849 121906 121957 121968 121981 122001 122058
[1366] 122103 122123 122132 122203 122204 122266 122297 122298 122299 122301 122302 122310 122327 122328 122329
[1381] 122330 122332 122375 122392 122415 122429 122439 122459 122460 122462 122464 122529 122567 122574 122580
[1396] 122581 122586 122605 122647 122655 122677 122681 122697 122718 122756 122758 122762 122763 122767 122768
[1411] 122780 122783 122784 122793 122808 122814 122815 122816 122858 122922 122925 122926 122927 122934 122937
[1426] 122938 122940 122953 122972 123032 123038 123039 123057 123059 123060 123061 123063 123064 123072 123085
[1441] 123103 123128 123210 123332 123350 123359 123403 123469 123519 123521 123522 123672 123675 123677 123709
[1456] 123739 123742 123762 123822 123923 123927 123952 123953 123980 123981 123994 124011 124014 124026 124075
[1471] 124078 124189 124215 124268 124785 124817 125046 125057 125197 125198 125200 125201 125203 125204 125205
[1486] 125395 125466 125468 125469 125474 125586 125619 125708 126067 126078 126110 126275 126276 126773 126774
[1501] 126799 126800 127171 127172 127611 127735 127780 127808 128105 128475 128837 128854 128868 128902 129014
[1516] 129523 129524 129980 130046 130059 130071 130074 130111 130118 130124 130129 130131 130134 130138 130142
[1531] 130149 130150 130151 130156 130160 130162 130165 130166 130167 130169 130170 130171 130172 130173 130174
[1546] 130175 130176 130178 130181 130183 130184 130187 130189 130190 130191 130192 130193 130194 130196 130207
[1561] 130208 130215 130219 130226 130227 130230 130231 130232 130238 130239 130241 130242 130247 130248 130249
[1576] 130252 130255 130261 130339 130347 130359 130420 130476 130763 130883 130886 130888 130889 130901 130902
[1591] 130903 131476 131480 131485 131910 131962 132293 132340 132457 132458 132618 132636 132638 132641 132642
[1606] 132889 133165 133217 133218 133401 133402 133407 133412 133415 133423 133532 133730 133733 133765 133772
[1621] 133773 133878 134308 134369 134377 134398 134847 134989 134994 135000 135001 135094 135275 135291 135294
[1636] 135305 135310 135314 135332 135335 135396 135397 135414 135415 135417 135418 135419 135423 135424 135452
[1651] 135455 135457 135458 135459 135462 135465 135489 135544 135938 135939 136097 136112 136201 136238 136240
[1666] 136281 136799 137915 137916 137941 138014 138081 138088 138093 138100 138169 138332 138533 138642 138643
[1681] 138645 138767 138768 138874 139205 139206 139442 139556 139580 140040 140150 140196 140209 140269 140279
[1696] 140553 140591 140595 140600 140601 140602 140604 140611 140614 140625 140626 140642 140655 140656 140670
[1711] 140697 140704 140810 140819 140825 140902 140962 141122 141291 141563 141564 141565 141576 141577 141585
[1726] 141603 141604 141792 141882 142086 142371 142400 142421 142428 142429 142511 142638 142777 142778 142779
[1741] 142780 142781 142797 142813 142833 142840 142927 143127 143184 143185 143220

$nb_outliers
[1] 1751

Outliers_Bayescan <- positive$SNP
length (Outliers_Bayescan)
#1236

######################################################## PCADAPT ###############################################

######### Plink VCF to bed to pcadapt ################


#1. Descargas plink y en la misma carpeta pones tu vcf y corres la siguiente linea:

VCF=/Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf

./plink2 --vcf $VCF --allow-extra-chr --make-bed --out /Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.bed

############## En R #######################

#Descargar la librerias

install.packages ("viridisLite")

library(adegenet)
library(pcadapt)
library(vcfR)
library(ggplot2)
library(dplyr)
library(viridis)

#Descargue los datos con la función read.pcadapt.

setwd("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/intSamples_withoutLD_1_100_K=8")

data <- read.pcadapt("intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.bed", type = "bed")

data2 <- read.vcfR("intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf")

#Scanning file to determine attributes.
#File attributes:
#  meta lines: 2190
#  header_line: 2191
#  variant count: 143272
#  column count: 109
#Meta line 2190 read in.
#All meta lines processed.
#gt matrix initialized.
#Character matrix gt created.
#  Character matrix gt rows: 143272
#  Character matrix gt cols: 109
#  skip: 0
#  nrows: 143272
#  row_num: 0
#Processed variant: 143272
#All variants processed
#2

data3 <- vcfR2genlight(data2)
snp <- as.data.frame(data3@loc.names)
ind <- as.data.frame(data3@ind.names)
colnames(ind) <-"INDIVIDUALS"

#Añadir información del mapa de población.

setwd("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/intSamples_withoutLD_1_100_K=8")

poplist.names <- c(rep("AIRE",1), rep("ESCURT", 1), rep("ADRG", 1), rep("MOLTONA", 1), rep("ESCURT",22), rep("MOLTONA",10), rep("GUARDIA",15), rep("AIRE",10), rep("COLOM",10), rep("PRCV",10), rep("ADRG",9), rep("ROVELLS",10))
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

#Para nuestro conjunto de datos, ejecutaremos pcadapt con K = 8 (ya que Se recomienda mantener los PC que corresponden a los valores propios a la izquierda de la línea recta (regla de Cattell).

#Por defecto, el parámetro min.maf está fijado en el 5%; aquí lo cambiamos al 1%. Los valores p de los SNPs con una frecuencia alélica menor que el umbral no se calculan (se devuelve NA).

data_pcadapt <- pcadapt(data, K = 8, min.maf = 0.01)

#Obtener una lista de valores P.
snps_pvalues <- cbind(snp, data_pcadapt$pvalues)
snps_pvalues_no_na <- na.omit(snps_pvalues)
write.table(snps_pvalues, "intsamples_withoutLD_Pvalues_K=8_.txt", sep="\t", quote=FALSE)

#Herramientas gráficas: el usuario también puede comprobar la distribución uniforme esperada de los valores p mediante un gráfico Q-Q

plot(data_pcadapt, option = "qqplot", col, snp.info = NULL, plt.pkg = "ggplot")

#Aquí el gráfico Q-Q confirma que la mayoría de los valores p siguen la distribución uniforme esperada.

#Compruebe la distribución de estos valores p.

hist(data_pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

#Un histograma de los valores p confirma que la mayoría de los valores p siguen una distribución uniforme. El exceso de valores P pequeños indica la presencia de valores atípicos.

#Visualizar la distribución de los valores p

quantile(snps_pvalues_no_na$`data_pcadapt$pvalues`, probs = c(0.01, 0.99))

#1%          99%
#1.763532e-11 9.999129e-01

#Obtenga sólo los marcadores que muestran valores p extremos: el 1% superior. #OJOOOOOOOO CAMBIAR EL NUMERO
#OJOOOOOOOOO

#top_1percent <- subset(snps_pvalues_no_na, snps_pvalues_no_na$`data_pcadapt$pvalues` <= 1.763532e-11)
#colnames(top_1percent) <- c("LOCUS","PVALUE")

#write.table(top_1percent, "intSamples_withoutLD_outliers_pcadapt_1_100_K=8.txt", sep="\t", quote=FALSE, row.names = FALSE)

##################### q-values to Cutoff  #######################

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("qvalue")

library("qvalue")
snps_pvalues_no_na$qval <- qvalue(snps_pvalues_no_na$`data_pcadapt$pvalues`)$qvalues
alpha <- 0.01
outliers_qvalue<- snps_pvalues_no_na$`data3@loc.names`[which(snps_pvalues_no_na$qval < alpha)]
length(outliers_qvalue)
#8892
write.table(outliers_qvalue, "intSamples_withoutLD_outliers_pcadapt_1_100_K=8.txt", sep="\t", quote=FALSE, row.names = FALSE)

####################################################### VENN DIAGRAM ###################################

#Cargar las librerias

#install.packages("VennDiagram")
library(VennDiagram)

#Guardar los nombres de SNP atípicos de Bayescan en un vector.
bayescan_outliers <- positive$SNP

#Compruebe el número de valores atípicos encontrados por BayeScan.

length(bayescan_outliers)
#[1] 1236

#Guarda los nombres de pcadapt en un vector.

pcadapt_outliers <- outliers_qvalue

#Comprueba el número de valores atípicos encontrados por pcadapt.

length(pcadapt_outliers)
#[1] 8892

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
write.table(outliers_comunes, "outliers_comunes_intSamples_1_100_K=8.txt", row.names=F, quote=F)
outliers_comunes

length(intersect(bayescan_outliers,pcadapt_outliers))
#1069

#Como se puede ver, el subconjunto del número compartido de candidatos es bajo, lo que sugiere que los métodos de exploración del genoma son muy variables.
#Este hallazgo está de acuerdo con Dalongeville et al. 2018, quienes mostraron que la identificación y el número de los valores atípicos son drásticamente diferentes dependiendo del método utilizado.

#Encuentre los valores atípicos únicos a BayeScan.

outliers_unicos_bayescan <- setdiff(bayescan_outliers,pcadapt_outliers)
write.table(outliers_unicos_bayescan, "outliers_unicos_bayescan_intSamples_1_100_K=8.txt", row.names=F, quote=F)
outliers_unicos_bayescan

#Ver el numero de Outliers enconrtrado por Bayescan

length(setdiff(bayescan_outliers,pcadapt_outliers))
#[1] 167

#Encuentre los valores atípicos únicos para pcadapta.

outliers_unicos_pcadapt <- setdiff(pcadapt_outliers,bayescan_outliers)
write.table(outliers_unicos_pcadapt, "outliers_unicos_pcadapt_intSamples_1_100_K=8.txt", row.names=F, quote=F)
outliers_unicos_pcadapt

#Ver el numero de Outliers enconrtrado por PcaAdapt
length(setdiff(pcadapt_outliers,bayescan_outliers))
#[1] 7823

#IR A VCFTOOLS Y EXTRAER LOS OUTLIERS DE BayeScan
# Hacer un archivo con los outrliers de bayescan con la informacion del cromosoma un TAB y luego la posicion :  outliers_comunes_intSamples_1_100_K=15_2.txt
#En vcftools la opcion positions me permite extraer los snps que me interesan

cd /Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/intSamples_withoutLD_1_100_K=8

vcftools --vcf intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf --positions outliers_comunes_intSamples_1_100_K=8_2.txt --recode --recode-INFO-all --out outliers_comunes_bayescan_pcadapt_intSamples_1_100_K=8

grep -v # outliers_comunes_bayescan_pcadapt_intSamples_1_100_K=8.recode.vcf | wc -l
