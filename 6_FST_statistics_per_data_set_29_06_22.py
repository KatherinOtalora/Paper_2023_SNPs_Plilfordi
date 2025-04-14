######################### FST ESTADISTICAS SOLO PARA LOS OUTLIERS COMUNES DE BAYESCAN Y PCADAPT###############
....................................ALLSAMPLES.......................................

# Vamos a Phyton
# En la terminal poner python y ya se entra automaticamente al entorno
#Entramos a la carpeta que tiene la informacion
import os
os.chdir(r'/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/allSamples_0.90_withoutLD_bayescan_1_100_K=15')

#Definimos nuestras variables
#OJO IMPORTANTE RECUERDA QUE EL ARCHIVO DE FST Y LOG PO ES UN OUTPUT QUE TE CALCULO BAYESCAN Y SE USARA PARA HACER ESTE FILTRO
# No olvidar poner SNP al nombre de la columna de los outliers comunes para que se pueda luego hacer el merge 
import pandas as pandas
df_a= pandas.read_csv("outliers_comunes_allSamples_0.90_1_100_K=15.txt")
df_b= pandas.read_csv("allSamples_0.90_filtred_withoutLD_bayescan_results_1_100_K=15.txt", delimiter='\t')

#Extraemos unicamente la informacion de los outliers comunes entre PCAdapt y Bayescan
result= pandas.merge(df_a, df_b, how='left')
print(result)

#Escribimos el nombre donde queremos que se guarde el output

writePath = 'outliers_comunes_allSamples_0.90_1_100_K=15_fst_statisc.txt'

with open(writePath, 'w') as f:
    df_string = result.to_string()
    f.write(df_string)
#Dar DOBLE enter al final para que ejecute la instruccion

....................................INTSAMPLES.......................................
# Vamos a Phyton
#Entramos a la carpeta que tiene la informacion
import os
os.chdir(r'/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/intSamples_withoutLD_1_100_K=8')

#Definimos nuestras variables
#OJO IMPORTANTE RECUERDA QUE EL ARCHIVO DE FST Y LOG PO ES UN OUTPUT QUE TE CALCULO BAYESCAN Y SE USARA PARA HACER ESTE FILTRO#OJO IMPORTANTE RECUERDA QUE EL ARCHIVO DE FST Y LOG PO ES UN OUTPUT QUE TE CALCULO BAYESCAN Y SE USARA PARA HACER ESTE FILTRO
# No olvidar poner SNP al nombre de la columna de los outliers comunes para que se pueda luego hacer el merge 
import pandas as pandas
df_a= pandas.read_csv("outliers_comunes_intSamples_1_100_K=8.txt")
df_b= pandas.read_csv("intSamples_0.90_filtred_withoutLD_bayescan_results_1_100_K=8.txt", delimiter='\t')

#Extraemos unicamente la informacion de los outliers comunes entre PCAdapt y Bayescan
result= pandas.merge(df_a, df_b, how='left')
print(result)

#Escribimos el nombre donde queremos que se guarde el output
writePath = 'outliers_comunes_intSamples_1_100_K=8_fst_statisc.txt'

with open(writePath, 'w') as f:
    df_string = result.to_string()
    f.write(df_string)
#Dar DOBLE enter al final para que ejecute la instruccion

....................................EXTSAMPLES.......................................
# Vamos a Phyton
#Entramos a la carpeta que tiene la informacion
import os
os.chdir(r'/Users/katherin_otalora/Documents/Bayescan_PCAdapt_Outliers/BayeScan2.1/binaries/extSamples_withoutLD_bayescan_1_100_K=10')

#Definimos nuestras variables
#OJO IMPORTANTE RECUERDA QUE EL ARCHIVO DE FST Y LOG PO ES UN OUTPUT QUE TE CALCULO BAYESCAN Y SE USARA PARA HACER ESTE FILTRO#OJO IMPORTANTE RECUERDA QUE EL ARCHIVO DE FST Y LOG PO ES UN OUTPUT QUE TE CALCULO BAYESCAN Y SE USARA PARA HACER ESTE FILTRO
# No olvidar poner SNP al nombre de la columna de los outliers comunes para que se pueda luego hacer el merge 
import pandas as pandas
df_a= pandas.read_csv("outliers_comunes_extSamples_1_100_K=10.txt")
df_b= pandas.read_csv("extSamples_0.90_filtred_withoutLD_bayescan_results_1_100_K=10.txt", delimiter='\t')

#Extraemos unicamente la informacion de los outliers comunes entre PCAdapt y Bayescan
result= pandas.merge(df_a, df_b, how='left')
print(result)

#Escribimos el nombre donde queremos que se guarde el output
writePath = 'outliers_comunes_extSamples_1_100_K=10_fst_statisc.txt'

with open(writePath, 'w') as f:
  df_string = result.to_string()
  f.write(df_string)
#Dar DOBLE enter al final para que ejecute la instruccion
