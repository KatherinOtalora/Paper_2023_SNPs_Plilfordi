##################### Outliers WITHOUTLD EXTSAMPLES K=8 #############################

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

cd /Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/extSamples_withoutLD_bayescan_1_100_K=10

#Corremos el Bayescan con 1:100
#Testeamos: 20 short pilots runs with 5000 integration, burn in set to 5 x 104 and thinning interval 10 and prior odds of 100.

/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/BayeScan2.1_macos64bits /Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/extSamples_withoutLD_bayescan_1_100_K\=10/extSamples_filtred_withoutLD_bayescan_1.txt  -od ./ -threads 8 \ -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100

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

setwd("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/extSamples_withoutLD_bayescan_1_100_K=10")

bayescan=read.table("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/extSamples_withoutLD_bayescan_1_100_K=10/extSamples_filtred_withoutLD_bayescan_1_fst.txt")

#La primera columna del marco de datos de bayescan es el ID del SNP.
#Las tres columnas siguientes (prob, log10(P0), y qval) están relacionadas con la prueba de adaptación local considerando el logaritmo de las probabilidades posteriores - log10(PO) - y el valor q para el modelo con selección.
#La quinta columna da el tamaño del efecto específico del locus (parámetro alfa). La última proporciona el FST específico del locus promediado en todas las poblaciones.

#Descargue la lista de SNPs en el orden correcto. El formato .geste tiene los SNPs en el mismo orden que el vcf utilizado para producir el formato .geste. Por lo tanto, puede utilizar este comando en bash para extraer la tercera columna que contiene la información de identificación de cada SNP en su vcf.
#Importante para el archivo orginal de vcf hay que cortar la columna 1 y 2 y luego juntarlas porque en la columna tres no hay informacion solo un . (Se juntan en algun editor de texto con _ entre medio de la culmna 1 y 2)

grep -v "#" extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf | cut -f 1,2 > extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_snps_1_100_K=10.txt

#A continuación, importe la información de los SNPs.

SNPb=read.table("extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_snps_1_100_K=10.txt", header=FALSE)

#Combinar los nombres de los valores atípicos con los resultados del marco de datos del bayescan.

bayescan=cbind(SNPb, bayescan)

#Cambiar el nombre de las columnas.

colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST")

#Escribe los resultados.

write.table(bayescan, "extSamples_0.90_filtred_withoutLD_bayescan_results_1_100_K=10.txt", quote=FALSE, sep="\t", row.names=FALSE)

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
#       84230         1606        81107

#Escriba los resultados de los SNPs potencialmente bajo selección (valor q < 0,01).

write.table(neutral, "extSamples_1_100_K=10_neutral.txt", row.names=F, quote=F)
write.table(balancing, "extSamples_1_100_K=10_balancing.txt", row.names=F, quote=F)
write.table(positive, "extSamples_1_100_K=10_positive.txt", row.names=F, quote=F)

#Transformación Log del valor Q para crear el gráfico ggplot.

range(bayescan$Q_VALUE)

#[1] 0.0001 0.9707

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

plot_bayescan("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/extSamples_withoutLD_bayescan_1_100_K=10/extSamples_filtred_withoutLD_bayescan_1_fst.txt")

$outliers
   [1]      3      4      5      7     11     14     17     18     20     28     30     31     34     36     49
  [16]     60     61    130    150    154    245    471    506    593    680    763    783    786    897    962
  [31]   1346   1937   1941   2028   2069   2203   2508   2529   2887   3279   3296   3717   3973   3974   4069
  [46]   4110   4111   4148   4288   4678   4690   4732   4782   4807   4828   4881   4944   5069   5072   5155
  [61]   5345   5432   5443   5487   5605   5629   5708   5979   5999   6025   6076   6087   6104   6105   6106
  [76]   6107   6109   6139   6223   6266   6270   6414   6415   6437   6455   6546   6549   6550   6748   6786
  [91]   6843   6882   6883   6884   6889   6896   7111   7278   7288   7321   7401   7402   7403   7404   7405
 [106]   7528   7530   7557   7649   7667   8140   8173   8174   8356   8447   8465   8479   8504   8505   8537
 [121]   8606   8608   8653   8701   8737   8959   8976   8987   8989   9073   9082   9098   9175   9288   9294
 [136]   9296   9500   9505   9518   9545   9560   9997  10065  10317  10368  10508  10606  10616  10641  10674
 [151]  10836  10877  10878  10924  10928  11127  11231  11479  11506  11643  11838  11994  12003  12022  12030
 [166]  12146  12463  12518  12588  12611  12620  12690  12694  12709  12710  12773  12774  12800  12804  12916
 [181]  13113  13177  13194  13223  13313  13315  13353  13354  13358  13510  13568  13680  13843  13949  14030
 [196]  14042  14047  14223  14291  14474  14746  14762  14763  14862  14868  15333  15335  15338  15339  15341
 [211]  15346  15348  15355  15358  15386  15389  15399  15400  15401  15403  15412  15413  15416  15417  15422
 [226]  15423  15426  15427  15445  15449  15484  15517  15566  15606  15610  15613  15617  15795  15852  15881
 [241]  15922  16045  16155  16189  16224  16330  16375  16413  16456  16471  16491  16499  16510  16699  16788
 [256]  16852  17009  17421  17615  17720  17970  18077  18078  18178  18624  18658  18762  18848  19020  19021
 [271]  19022  19024  19025  19039  19046  19155  19195  19273  19275  19298  19299  19300  19334  19444  19659
 [286]  19660  20016  20027  20152  20182  20247  20248  20378  20457  20540  20566  20617  20634  20640  20656
 [301]  20725  20727  20825  20826  20866  20947  20992  21081  21116  21162  21166  21183  21197  21243  21292
 [316]  21309  21343  21445  21453  21534  21537  21538  21539  21611  22041  22068  22159  22160  22161  22162
 [331]  22163  22164  22165  22166  22544  22607  22694  22793  22928  23125  23250  23255  23282  23284  23339
 [346]  23340  23594  23627  23629  23631  23749  23831  23878  23922  23957  23974  23987  24112  24164  24166
 [361]  24195  24226  24297  24336  24339  24411  24437  24850  24996  25006  25067  25100  25151  25164  25207
 [376]  25211  25214  25220  25256  25406  25418  25632  25879  25943  25948  26042  26290  26305  26306  26425
 [391]  26514  26753  26803  26829  27268  27733  27736  27848  27948  28100  28671  28769  28862  28872  28873
 [406]  28919  28922  28931  28934  28943  28947  28953  28958  28976  28991  28994  29071  29116  29126  29133
 [421]  29134  29135  29153  29161  29164  29165  29187  29350  29585  29668  29841  29930  29931  29932  29938
 [436]  29939  29941  29975  29989  30344  30376  30377  30815  31066  31252  31438  31567  31747  32271  32327
 [451]  32328  32538  32649  32675  32885  32886  32917  33070  33219  33340  33402  33562  33565  33587  33735
 [466]  33967  34257  34308  34574  34612  34635  34850  34900  35140  35141  35142  35143  35144  35263  35280
 [481]  35380  35887  36185  36237  36321  36324  36352  36357  36393  36546  36551  36573  36577  36580  36866
 [496]  36928  36984  36985  36987  37485  37489  37537  37611  37658  37757  37791  37971  38045  38094  38171
 [511]  38175  38301  38387  38406  38410  38443  38457  38465  38519  38566  38576  38592  38682  38784  38831
 [526]  39338  39366  39462  39491  39492  39494  39809  39843  39857  40376  40395  40420  40613  40614  40616
 [541]  40663  40683  40727  40734  40762  40820  40839  41170  41368  41369  41782  41994  42036  42334  42338
 [556]  42339  42426  42427  42476  42477  42712  42769  42778  42901  42968  42998  43048  43140  43429  43453
 [571]  43481  43532  43587  43600  43661  43752  43900  44407  45029  45114  45413  45414  45417  45896  46078
 [586]  46122  46315  46316  46317  46453  46481  46598  46599  46627  46639  46807  47094  47095  47098  47281
 [601]  47299  47336  47378  47599  47699  47704  47727  47747  47751  47759  47802  47855  47912  47953  48069
 [616]  48246  48247  48321  48463  48494  48639  49079  49331  49426  49427  49589  49797  49798  49929  49930
 [631]  49986  50011  50017  50064  50135  50287  50302  50303  50304  50306  50307  50318  50319  50406  50473
 [646]  50523  50736  50815  50867  51068  51210  51286  51369  51521  51736  51746  51797  51975  52000  52382
 [661]  52478  52481  52580  52857  53164  53409  53421  53425  53691  53692  54079  54106  54107  54439  54444
 [676]  54468  54669  54705  54707  54724  54782  54785  54786  54787  54789  54790  54808  54812  54813  54815
 [691]  54819  54823  54825  54836  54837  54838  54843  54845  54848  54859  54861  54870  54897  54912  54913
 [706]  55012  55208  55312  55515  55735  55736  55779  55793  55804  55860  56050  56143  56358  56390  56512
 [721]  56515  56523  56635  56654  56899  56928  57058  57320  57394  57395  57428  57480  57560  57822  57935
 [736]  57936  58156  58477  58649  58703  58998  59054  59138  59289  59508  59575  59586  59599  59694  59739
 [751]  59975  59976  60038  60039  60084  60127  60378  60379  60394  60459  60462  60549  60570  60713  60737
 [766]  61011  61015  61090  61091  61092  61173  61304  61517  61575  61577  61578  61579  61580  61878  61897
 [781]  61933  62012  62028  62038  62039  62231  62233  62275  62327  62455  62593  62728  62733  62833  62958
 [796]  62964  63093  63097  63106  63107  63173  63206  63244  63350  63352  63353  63356  63838  64029  64197
 [811]  64407  64409  64443  64489  64608  64836  64954  64970  65002  65011  65252  65255  65256  65332  65361
 [826]  65362  65364  65365  65366  65553  65584  65601  65678  65864  65951  65952  65953  65977  65980  65987
 [841]  66002  66006  66024  66042  66056  66061  66067  66068  66073  66074  66075  66077  66085  66361  66362
 [856]  66533  66535  66538  66587  66599  66694  67015  67041  67240  67463  67798  67899  67993  68048  68212
 [871]  68412  68518  68592  68678  68781  68841  69105  69252  69324  69350  69410  69431  69741  69927  69930
 [886]  70181  70182  70330  70339  70511  70556  70612  70673  70789  70792  70833  71090  71186  71313  71314
 [901]  71316  71318  71321  71646  71909  72302  72338  72489  72816  72817  72821  72879  72880  73118  73176
 [916]  73177  73465  73514  73564  73694  73706  74059  74067  74105  74291  74318  74325  74604  74619  74740
 [931]  75162  75205  75283  75468  75472  75584  75687  75728  76005  76008  76123  76124  76127  76129  76130
 [946]  76172  76334  76394  76449  76486  76504  76505  76626  77002  77146  77168  77171  77275  77311  77457
 [961]  77529  77568  77577  77655  77824  78105  78182  78251  78353  78473  78739  78821  78977  79317  79359
 [976]  79595  79797  79799  79800  79854  79871  79997  80105  80133  80204  80308  80316  80486  80545  80648
 [991]  80652  80653  80716  80734  80928  80952  81133  81137  81161  81247  81353  81410  81422  81494  81600
[1006]  81614  81850  81851  81852  81853  81854  81856  81881  81985  82054  82102  82117  82242  82399  82422
[1021]  82423  82743  82910  82911  82932  82933  83062  83319  83320  83425  83490  83530  84003  84044  84084
[1036]  84201  84593  84594  84645  84776  84855  85000  85042  85304  85309  85412  85413  85736  85796  85959
[1051]  85961  85962  86048  86252  86330  86373  86502  86504  86529  86998  87063  87087  87142  87240  87273
[1066]  87289  87300  87301  87302  87318  87329  87334  87337  87338  87339  87342  87347  87349  87350  87353
[1081]  87356  87357  87367  87374  87376  87388  87393  87395  87396  87397  87398  87454  87491  87717  87720
[1096]  87721  87858  87947  88636  88872  88939  88941  88958  88960  89090  89299  89348  89349  89350  89351
[1111]  89571  90040  90041  90186  90292  90296  90301  90302  90306  90311  90359  90388  90389  90460  90507
[1126]  90508  90546  90581  90674  90675  90781  90911  90982  91157  91166  91266  91461  91534  91579  92227
[1141]  92233  92234  92387  92388  92709  92722  92773  92880  92948  93066  93089  93097  93151  93153  93198
[1156]  93199  93200  93203  93284  93286  93407  93773  93774  93800  93968  94223  94226  94233  94304  94305
[1171]  94306  94380  94381  94382  94624  94652  94710  94787  94825  94881  94908  95039  95114  95115  95312
[1186]  95352  95482  95674  95918  95923  95940  95959  95961  95962  95968  96028  96029  96030  96031  96036
[1201]  96742  96852  97373  97379  97382  97384  97391  97406  97407  97408  97411  97443  97511  97721  97949
[1216]  97954  97958  97962  97974  97979  97982  97983  97989  97998  98000  98001  98020  98032  98045  98578
[1231]  98869  98893  98906  98975  99043  99180  99182  99348  99701  99722  99723  99743  99879  99882 100044
[1246] 100062 100738 100788 100795 100814 101105 101107 101129 101295 101323 101430 101528 101552 101860 101987
[1261] 102033 102035 102169 102178 102188 102292 102313 102380 102381 102584 102606 102665 102797 102984 103199
[1276] 103216 103255 103319 103370 103415 103841 103851 103875 103936 104190 104337 104348 104349 104371 104378
[1291] 104461 104571 104711 104866 105437 105713 105725 105830 106196 106224 106225 106226 106231 106272 106317
[1306] 106442 106462 106470 106493 106500 106517 106562 106669 106670 106671 106672 106673 106675 106676 106677
[1321] 106680 106681 106685 106691 106706 106708 106721 106723 106727 106744 106749 106759 106760 106776 106777
[1336] 106856 106871 106872 106961 106987 107072 107233 107411 108011 108171 108172 108186 108284 108359 108484
[1351] 108576 108614 108736 108738 108830 108845 109068 109077 109398 109514 109521 109701 109756 110240 110242
[1366] 110475 110502 110532 110548 110566 110640 110642 110797 110838 110991 110992 110993 110998 111044 111046
[1381] 111047 111354 111423 111666 111670 111795 111797 112056 112058 112106 112121 112122 112217 112487 112565
[1396] 112668 112669 112743 112892 112895 112951 112953 112983 113106 113148 113291 113364 113500 113501 114167
[1411] 114177 114190 114207 114212 114266 114518 114603 114644 114646 114647 114649 114652 114672 114674 114679
[1426] 114680 114681 114682 114686 114688 114689 114690 114691 114700 114701 114711 114712 114717 114727 114740
[1441] 114742 114792 114796 114800 114813 114826 114906 114914 114916 114917 114961 114965 115110 115337 115351
[1456] 115409 115570 115587 115806 116243 116294 116295 116409 116508 116509 116692 116695 116875 116951 117543
[1471] 117544 117545 117546 118045 118126 118127 118230 118525 118540 118543 118740 118741 118814 118833 118834
[1486] 118871 118932 118982 118993 119103 119113 119121 119148 119273 119322 119349 119540 119583 119842 119866
[1501] 120012 120446 120513 120517 120523 120528 120629 120762 120911 121054 121115 121154 121215 121233 121259
[1516] 121329 121560 121667 121727 121778 121839 121857 121858 121859 121861 121862 121863 121864 121865 121867
[1531] 121868 121870 121873 121888 121980 122212 122251 122433 122512 122624 122719 122835 122836 122867 122868
[1546] 122955 123219 123248 123277 123449 123452 123545 123721 123914 124210 124328 124391 124466 124487 124978
[1561] 124982 124983 124991 124992 125060 125160 125249 125511 125516 125517 125685 125843 126048 126142 126236
[1576] 126270 126353 126369 126418 126478 126556 126693 126777 126779 126821 127051 127159 127272 127336 128053
[1591] 128231 128336 128356 128357 128380 128382 128389 128406 128415 128421 128423 128424 128425 128429 128431
[1606] 128448 128451 128452 128461 128464 128473 128478 128482 128484 128487 128494 128495 128501 128801 128869
[1621] 128892 129106 129107 129149 129408 129482 129490 129564 129691 129767 129792 129862 129906 129910 130005
[1636] 130318 130427 130841 130849 130850 130851 130891 131330 131379 131460 131461 131483 131613 131690 131691
[1651] 131850 131897 132048 132056 132096 132105 132106 132107 132108 132110 132201 132359 132654 132822 132895
[1666] 133001 133040 133072 133152 133453 133632 133667 133699 133868 133891 134022 134023 134319 134456 134570
[1681] 134578 134828 134829 134965 134970 135219 135231 135435 135470 136108 136116 136124 136168 136246 136388
[1696] 136592 136820 136965 137082 137099 137116 137142 137179 137293 137295 137321 137578 137739 137871 138059
[1711] 138162 138234 138269 138295 138298 138381 138382 138685 138816 138817 138818 138821 138822 138843 138851
[1726] 138852 138983 139002 139081 139082 139187 139314 139315 139425 139482 139560 139566 139665 140476 140478
[1741] 140480 140489 140491 141060 141151 141199 141267 141306 141307 141635 141778 141850 141851 141852 141879
[1756] 141880 141904 141972 141977 142015 142111 142153 142154 142155 142156 142161 142164 142174 142210 142243
[1771] 142253 142254 142262 142271 142281 142288 142289 142290 142292 142293 142294 142295 142296 142297 142298
[1786] 142299 142300 142301 142302 142303 142304 142306 142307 142308 142309 142310 142313 142314 142315 142316
[1801] 142317 142318 142319 142320 142321 142322 142324 142325 142326 142327 142328 142329 142331 142332 142334
[1816] 142337 142339 142347 142351 142356 142357 142358 142367 142372 142374 142375 142376 142381 142382 142384
[1831] 142385 142389 142392 142393 142394 142400 142402 142405 142407 142410 142414 142419 142420 142421 142426
[1846] 142427 142428 142429 142430 142433 142434 142435 142436 142442 142443 142452 142462 142465 142466 142473
[1861] 142488 142491 142494 142495 142497 142498 142506 142507 142509 142510 142515 142521 142524 142527 142528
[1876] 142535 142536 142543 142549 142550 142573 142628 142636 142770 142779 142782 142786 142788 142826 142843
[1891] 142844 142861 142891 142895 142971 142981 142982 142992 143003 143006 143007 143009 143011 143014 143015
[1906] 143022 143040 143043 143045 143065 143078 143096 143115 143138 143156 143158 143159 143186 143201 143231
[1921] 143240 143242 143244 143268 143293 143297 143308 143348 143363 143382 143478 143504 143529 143592 143671
[1936] 143672 143679 143730 143778 143793 143799 143833 143861 143875 143894 143895 143896 143898 143902 143903
[1951] 143913 143944 143946 144029 144134 144135 144140 144141 144147 144150 144153 144156 144160 144161 144189
[1966] 144190 144222 144223 144224 144247 144267 144270 144276 144279 144281 144300 144302 144362 144367 144368
[1981] 144369 144427 144444 144454 144497 144505 144543 144544 144567 144571 144602 144621 144637 144651 144655
[1996] 144675 144676 144698 144730 144733 144743 144744 144763 144766 144790 144804 144812 144841 144847 144848
[2011] 144869 144872 144873 144885 144935 144980 145005 145128 145133 145216 145249 145275 145277 145285 145321
[2026] 145324 145399 145400 145401 145402 145461 145488 145489 145494 145495 145496 145538 145542 145562 145567
[2041] 145575 145590 145669 145715 145717 145732 145733 145811 145821 145831 145832 145835 145858 145965 146116
[2056] 146125 146126 146149 146302 146337 146377 146431 146515 146582 146593 146594 146619 146687 146688 146745
[2071] 146768 146787 146805 146887 146927 146990 147039 147112 147117 147119 147120 147121 147122 147123 147167
[2086] 147168 147173 147174 147236 147252 147273 147274 147277 147287 147294 147295 147296 147297 147298 147299
[2101] 147373 147463 147468 147469 147534 147690 147693 147701 147711 147721 147829 147882 147883 147937 147939
[2116] 147947 147948 147949 147950 147998 148229 148333 148348 148597 148600 148734 148817 148876 148944 149035
[2131] 149075 149242 149341 149711 149742 150056 150135 150271 150272 150514 150543 150574 150614 150738 150766
[2146] 151132 151175 151190 151206 151207 151230 151423 151426 151497 151568 151569 151625 151634 151799 151800
[2161] 151827 151920 152397 152428 152450 152452 152601 152690 152756 152757 152769 152819 152842 152844 152857
[2176] 152858 152860 152869 152888 152895 152910 152913 152914 152920 152925 152926 152936 152939 152940 152943
[2191] 152944 152945 152946 152949 152952 152953 152954 152955 152957 152959 152960 152961 152963 152964 152965
[2206] 152966 152967 152968 152969 152970 152971 152973 152974 152992 152999 153000 153001 153002 153003 153005
[2221] 153012 153015 153031 153056 153057 153058 153059 153060 153062 153067 153084 153148 153159 153264 153299
[2236] 153301 153403 153404 153519 153541 153622 153718 153722 153781 154030 154031 154080 154081 154220 154503
[2251] 154536 154537 154683 154684 154685 154692 154869 154942 154968 155031 155210 155211 155252 155265 155319
[2266] 155533 155609 155979 156093 156479 157150 157338 157339 157340 157378 157487 157488 157524 157934 157935
[2281] 157944 158305 158339 158353 158457 158529 158781 158782 158783 158784 158785 158913 158921 158922 158924
[2296] 158944 158946 158947 158964 158969 158971 158976 158995 159003 159016 159019 159025 159036 159040 159050
[2311] 159057 159060 159061 159065 159093 159159 159230 159231 159232 159234 159235 159236 159416 159494 159607
[2326] 159726 159730 159731 159742 159743 159744 159745 159754 159759 159832 159855 160392 160443 160444 160445
[2341] 160446 160463 160464 160469 160474 160484 160491 160524 160532 160533 160535 160536 160537 160555 160557
[2356] 160558 160559 160561 160623 160798 160889 160938 161367 161442 161477 161507 161508 161647 161726 161833
[2371] 162879 162981 163000 163032 163321 163363 163460 163593 163660 163828 163865 163930 164054 164063 164270
[2386] 164271 164282 164302 164318 164358 164374 164619 164620 164871 164994 165099 165172 165456 165469 165483
[2401] 165485 165489 165491 165534 165583 165718 165804 165805 165806 165826 165893 165918 165944 166028 166030
[2416] 166124 166128 166146 166211 166230 166240 166265 166270 166273 166274 166291 166326 166350 166360 166461
[2431] 166512 166559 166729 166816 166828 166872

$nb_outliers
[1] 2436

Outliers_Bayescan <- positive$SNP
length (Outliers_Bayescan)
#1606

######################################################## PCADAPT ###############################################

########## Plink VCF to bed to pcadapt ################


#1. Descargas plink y en la misma carpeta pones tu vcf y corres la siguiente linea:

VCF=/Users/katherin_otalora/Documents/Plink_1/plink/extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf

./plink2 --vcf $VCF --allow-extra-chr --make-bed --out /Users/katherin_otalora/Documents/Plink_1/plink/extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.bed

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

setwd("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/extSamples_withoutLD_bayescan_1_100_K=10")

data <- read.pcadapt("extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.bed", type = "bed")

data2 <- read.vcfR("extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf")

#Scanning file to determine attributes.
#File attributes:
#  meta lines: 2190
#  header_line: 2191
#  variant count: 166943
#  column count: 100
#Meta line 2190 read in.
#All meta lines processed.
#gt matrix initialized.
#Processed variant: 166943
#All variants processed
#3

data3 <- vcfR2genlight(data2)
snp <- as.data.frame(data3@loc.names)
ind <- as.data.frame(data3@ind.names)
colnames(ind) <-"INDIVIDUALS"

#Añadir información del mapa de población.

setwd("/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/extSamples_withoutLD_bayescan_1_100_K=10")

poplist.names <- c(rep("REI",1), rep("AIRE", 1), rep("REI", 9), rep("PORROS", 1), rep("AIRE",1), rep("PORROS",10), rep("AIRE",1), rep("FORADADA",10), rep("AIRE",1), rep("ESCLTS",10), rep("AIRE",1), rep("DRAGONERA",4), rep("CABRERAH",6), rep("AIRE",1), rep("CABRERAH",6), rep("COLOM",4), rep("AIRE",1), rep("COLOM",6), rep("CABREBRAI",4), rep("AIRE",1), rep("CABRERAI",4), rep("COLOMER",4), rep("AIRE",4))
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

#Para nuestro conjunto de datos, ejecutaremos pcadapt con K = 10 (ya que Se recomienda mantener los PC que corresponden a los valores propios a la izquierda de la línea recta (regla de Cattell).

#Por defecto, el parámetro min.maf está fijado en el 5%; aquí lo cambiamos al 1%. Los valores p de los SNPs con una frecuencia alélica menor que el umbral no se calculan (se devuelve NA).

data_pcadapt <- pcadapt(data, K = 10, min.maf = 0.01)

#Obtener una lista de valores P.
snps_pvalues <- cbind(snp, data_pcadapt$pvalues)
snps_pvalues_no_na <- na.omit(snps_pvalues)
write.table(snps_pvalues, "extsamples_withoutLD_Pvalues_K=10_.txt", sep="\t", quote=FALSE)

#Herramientas gráficas: el usuario también puede comprobar la distribución uniforme esperada de los valores p mediante un gráfico Q-Q

plot(data_pcadapt, option = "qqplot", col, snp.info = NULL, plt.pkg = "ggplot")

#Aquí el gráfico Q-Q confirma que la mayoría de los valores p siguen la distribución uniforme esperada.

#Compruebe la distribución de estos valores p.

hist(data_pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

#Un histograma de los valores p confirma que la mayoría de los valores p siguen una distribución uniforme. El exceso de valores P pequeños indica la presencia de valores atípicos.

#Visualizar la distribución de los valores p

quantile(snps_pvalues_no_na$`data_pcadapt$pvalues`, probs = c(0.01, 0.99))

#1%          99%
#1.516957e-13 9.986301e-01

#Obtenga sólo los marcadores que muestran valores p extremos: el 1% superior. #OJOOOOOOOO CAMBIAR EL NUMERO
#OJOOOOOOOOO

#top_1percent <- subset(snps_pvalues_no_na, snps_pvalues_no_na$`data_pcadapt$pvalues` <= 1.516957e-13)
#colnames(top_1percent) <- c("LOCUS","PVALUE")

#write.table(top_1percent, "extSamples_withoutLD_outliers_pcadapt_1_100_K=10.txt", sep="\t", quote=FALSE, row.names = FALSE)

##################### q-values to Cutoff  #######################

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("qvalue")

library("qvalue")
snps_pvalues_no_na$qval <- qvalue(snps_pvalues_no_na$`data_pcadapt$pvalues`)$qvalues
alpha <- 0.01
outliers_qvalue<- snps_pvalues_no_na$`data3@loc.names`[which(snps_pvalues_no_na$qval < alpha)]
length(outliers_qvalue)
#12740
write.table(outliers_qvalue, "extSamples_withoutLD_outliers_pcadapt_1_100_K=10.txt", sep="\t", quote=FALSE, row.names = FALSE)


####################################################### VENN DIAGRAM ###################################

#Cargar las librerias

install.packages("VennDiagram")
library(VennDiagram)

#Guardar los nombres de SNP atípicos de Bayescan en un vector.
bayescan_outliers <- positive$SNP

#Compruebe el número de valores atípicos encontrados por BayeScan.

length(bayescan_outliers)
#[1] 1606

#Guarda los nombres de pcadapt en un vector.

pcadapt_outliers <- outliers_qvalue

#Comprueba el número de valores atípicos encontrados por pcadapt.

length(pcadapt_outliers)
#[1] 12740

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
write.table(outliers_comunes, "outliers_comunes_extSamples_1_100_K=10.txt", row.names=F, quote=F)
outliers_comunes
length(intersect(bayescan_outliers,pcadapt_outliers))
#1271

#Como se puede ver, el subconjunto del número compartido de candidatos es bajo, lo que sugiere que los métodos de exploración del genoma son muy variables.
#Este hallazgo está de acuerdo con Dalongeville et al. 2018, quienes mostraron que la identificación y el número de los valores atípicos son drásticamente diferentes dependiendo del método utilizado.

#Encuentre los valores atípicos únicos a BayeScan.

outliers_unicos_bayescan <- setdiff(bayescan_outliers,pcadapt_outliers)
write.table(outliers_unicos_bayescan, "outliers_unicos_bayescan_extSamples_1_100_K=10.txt", row.names=F, quote=F)
outliers_unicos_bayescan

#Ver el numero de Outliers enconrtrado por Bayescan

length(setdiff(bayescan_outliers,pcadapt_outliers))
#[1] 335

#Encuentre los valores atípicos únicos para pcadaptarse.

outliers_unicos_pcadapt <- setdiff(pcadapt_outliers,bayescan_outliers)
write.table(outliers_unicos_pcadapt, "outliers_unicos_pcadapt_extSamples_1_100_K=10.txt", row.names=F, quote=F)
outliers_unicos_pcadapt

#Ver el numero de Outliers enconrtrado por Bayescan
length(setdiff(pcadapt_outliers,bayescan_outliers))
#[1] 11469

#IR A VCFTOOLS Y EXTRAER LOS OUTLIERS DE BayeScan
# Hacer un archivo con los outrliers de bayescan con la informacion del cromosoma un TAB y luego la posicion :  outliers_comunes_allSamples_0.90_1_100_K=15_2.txt
#En vcftools la opcion positions me permite extraer los snps que me interesan

cd /Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/extSamples_withoutLD_bayescan_1_100_K=10

vcftools --vcf extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf --positions outliers_comunes_extSamples_1_100_K=10_2.txt --recode --recode-INFO-all --out outliers_comunes_bayescan_pcadapt_extSamples_1_100_K=10

grep -v # outliers_comunes_bayescan_pcadapt_extSamples_1_100_K=10.recode.vcf | wc -l
