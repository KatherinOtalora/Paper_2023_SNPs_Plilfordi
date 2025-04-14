#..................................................................................................................................................................................................................................................................
#................................###### IMEDEA PROJECT AND BASSITA ET AL ###### .INTSAMPLES, EXTSAMPLES AND ALLSAMPLES.................................................................................................................................................................................................................................
#..................................................................................................................................................................................................................................................................

### NON FILTERED VCF FILES: COUNT NUMBER OF SNPS

bcftools #Para saber si esta correctamente instalado el programa

cd /Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/INT_SAMPLES

bcftools view -H /Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/INT_SAMPLES/intSample.g.geno.vcf.gz | wc -l #1.888.392

cd /Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/EXT_SAMPLES

bcftools view -H /Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/INT_SAMPLES/extSample.g.geno.vcf.gz  | wc -l #4.851.070

cd /Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/ALL_SAMPLES

bcftools view -H /Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/INT_SAMPLES/allSample.g.geno.vcf.gz  | wc -l #6.394.354

### FILTERING SNPS

#FILTERS (estos son comunes, no hay que repetirlos)

#..................................................................................................................................................................................................................................................................
#............................................................................INTSAMPLES......................................................................................................................................................................................
#..................................................................................................................................................................................................................................................................

VCF_IN=/Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/INT_SAMPLES/intSample.g.geno.vcf.gz
VCF_OUT=/Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/INT_SAMPLES/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf.gz


MINALLELE=2
MAXALLELE=2
MAF=0.05
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50
MAX_MISSING=0.75  #Este es un filtro que va a variar dependiendo de los analysis que queremos realizar

vcftools --gzvcf $VCF_IN --remove-indels --min-alleles $MINALLELE --max-alleles $MAXALLELE --maf $MAF --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --max-missing $MAX_MISSING --recode --stdout | gzip -c > $VCF_OUT

cat out.log

bcftools view -H intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf.gz | wc -l



#..................................................................................................................................................................................................................................................................
#................................................................................EXTSAMPLES..................................................................................................................................................................................
#..................................................................................................................................................................................................................................................................


VCF_IN=/Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/INT_SAMPLES/extSample.g.geno.vcf.gz
VCF_OUT=/Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/INT_SAMPLES/extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf.gz


MINALLELE=2
MAXALLELE=2
MAF=0.05
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50
MAX_MISSING=0.75  #Este es un filtro que va a variar dependiendo de los analisis que queremos realizar

vcftools --gzvcf $VCF_IN --remove-indels --min-alleles $MINALLELE --max-alleles $MAXALLELE --maf $MAF --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --max-missing $MAX_MISSING --recode --stdout | gzip -c > $VCF_OUT

cat out.log

bcftools view -H extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf.gz | wc -l

#..................................................................................................................................................................................................................................................................
#................................................................................ALLSAMPLES..MAX_MISSING=0.95................................................................................................................................................................................
#..................................................................................................................................................................................................................................................................


VCF_IN=/Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/ALL_SAMPLES/allSample.g.geno.vcf
VCF_OUT=/Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/ALL_SAMPLES/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095.vcf.gz



MINALLELE=2
MAXALLELE=2
MAF=0.05
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50
MAX_MISSING=0.95  #Este es un filtro que va a variar dependiendo de los analisis que queremos realizar
#MAX_MISSING=0.90 Se realizo tambien con esta opcion

vcftools --gzvcf $VCF_IN --remove-indels --min-alleles $MINALLELE --max-alleles $MAXALLELE --maf $MAF --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --max-missing $MAX_MISSING --recode --stdout | gzip -c > $VCF_OUT

cat out.log

bcftools view -H allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095.vcf.gz | wc -l

#Outputting VCF file...
#After filtering, kept 3418 out of a possible 6394354 Sites

#Filtro por Missindv

vcftools --gzvcf /Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/ALL_SAMPLES/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095.vcf.gz  --missing-indv --out allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_MISSINDV

#OUTPUT:".imiss"

#SACAR LOS INDIVIDUOS CON With more than 90% or 70%, 10% of missing data, #Para deshacerse de los individuos que no se secuenciaron bien, con mas del 90% de missing data

awk '$5 > 0.9' allSample.g.geno_FILTERED_MISSINDV.imiss| cut -f1 > lowDP_0.9.indv
awk '$5 > 0.7' allSample.g.geno_FILTERED_MISSINDV.imiss| cut -f1 > lowDP_0.7.indv
awk '$5 > 0.5' allSample.g.geno_FILTERED_MISSINDV.imiss| cut -f1 > lowDP_0.5.indv
awk '$5 > 0.4' allSample.g.geno_FILTERED_MISSINDV.imiss| cut -f1 > lowDP_0.4.indv

#..................................................................................................................................................................................................................................................................
#............................................................................ALL SAMPLES MAX_MISSING=0.90......................................................................................................................................................................................


VCF_IN=/Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/ALL_SAMPLES/allSample.g.geno.vcf
VCF_OUT=/Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/ALL_SAMPLES/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.vcf.gz



MINALLELE=2
MAXALLELE=2
MAF=0.05
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50
MAX_MISSING=0.90  #Este es un filtro que va a variar dependiendo de los analisis que queremos realizar

vcftools --gzvcf $VCF_IN --remove-indels --min-alleles $MINALLELE --max-alleles $MAXALLELE --maf $MAF --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --max-missing $MAX_MISSING --recode --stdout | gzip -c > $VCF_OUT

cat out.log

bcftools view -H /Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/Data_set_2022/ALL_SAMPLES/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.vcf.gz | wc -l

#Outputting VCF file...
#After filtering, kept 5523 out of a possible 6394354 Sites
#Run Time = 323.00 seconds

#..................................................................................................................................................................................................................................................................
#..................................................................................................................................................................................................................................................................
#.....................................................................................#### LINKAGE DISEQUILIBRIUM PRUNING #### ..........................................................................................................................................................................................
#..................................................................................................................................................................................................................................................................

#..................................................................................................................................................................................................................................................................
#...................................................................................########PLINK INTSAMPLES intSamples #####...............................................................................................................................................................................
#.........................................................................................................50_5_02.........................................................................................................................................................

cd /Users/katherin_otalora/Documents/Plink_1/plink

VCF=/Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf

#1. Descargas plink y en la misma carpeta pones tu vcf y corres la siguiente linea:

./plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 5 0.20 --out intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_prune_50_5_02

#2. PODAR EL ARCHIVO DE LOS SNPS CON LD:

./plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --extract intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_prune_50_5_02.prune.in --make-bed --pca --out /Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02

#3.TRANSFORMAR EL PRUNEDDATA ( ARCHIVO PODADO) A VCF DE NUEVO

./plink2  --bfile /Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02 --double-id --allow-extra-chr --set-missing-var-ids @:# --cow --recode vcf --out /Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02

grep -v "^#" /Users/katherin_otalora/Documents/Plink_1/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.vcf | wc -l


#..................................................................................................................................................................................................................................................................
#...................................................................................########PLINK EXTSAMPLES extSamples #####...............................................................................................................................................................................
#.......................................................................................................50_5_02...........................................................................................................................................................

cd /Users/katherin_otalora/Documents/Plink_1/plink

VCF=/Users/katherin_otalora/Documents/Plink_1/plink/extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf

#1. Descargas plink y en la misma carpeta pones tu vcf y corres la siguiente linea:

./plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 5 0.20 --out extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_prune_50_5_02

#2. PODAR EL ARCHIVO DE LOS SNPS CON LD:

./plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --extract extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_prune_50_5_02.prune.in --make-bed --pca --out /Users/katherin_otalora/Documents/Plink_1/plink/extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02

#3.TRANSFORMAR EL PRUNEDDATA ( ARCHIVO PODADO) A VCF DE NUEVO

./plink2  --bfile /Users/katherin_otalora/Documents/Plink_1/plink/extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02 --double-id --allow-extra-chr --set-missing-var-ids @:# --cow --recode vcf --out /Users/katherin_otalora/Documents/Plink_1/plink/extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02

grep -v "^#" /Users/katherin_otalora/Documents/Plink_1/plink/extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.vcf | wc -l

#Nos quedamos finalmente con el output de 50 5 0.20 con un total de 58207 SNPs para extSample y	34065 para intSample

#transform to .gz

bcftools view -O z -o extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.vcf.gz extSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075_LD_prunned_50_5_02.vcf | wc -l


#..................................................................................................................................................................................................................................................................
#...................................................................................########PLINK ALLSAMPLES allSamples maxMISS0.90 #####...............................................................................................................................................................................
#.......................................................................................................50_5_02...........................................................................................................................................................

VCF=/Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.vcf

cd /Users/katherin_otalora/Documents/Plink_1/plink

grep -v "^#" /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090.vcf | wc -l

#1. Descargas plink y en la misma carpeta pones tu vcf y corres la siguiente linea:

./plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 5 0.20 --out allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_prune_50_5_02

#2. PODAR EL ARCHIVO DE LOS SNPS CON LD:

./plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --extract allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_prune_50_5_02.prune.in --make-bed --pca --out /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02

#3.TRANSFORMAR EL PRUNEDDATA ( ARCHIVO PODADO) A VCF DE NUEVO

./plink2  --bfile /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02 --double-id --allow-extra-chr --set-missing-var-ids @:# --cow --recode vcf --out /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02

grep -v "^#" /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.vcf | wc -l


#Nos quedamos finalmente con el output de 50 5 0.20 con un total de 2876 SNPs para allSamples maxMISS0.90

#transform to .gz

bcftools view -O z -o allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.vcf.gz allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.vcf | wc -l


#..................................................................................................................................................................................................................................................................
#...................................................................................########PLINK ALLSAMPLES allSamples maxMISS0.95 #####...............................................................................................................................................................................
#.......................................................................................................50_5_02...........................................................................................................................................................

VCF=/Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095.vcf

cd /Users/katherin_otalora/Documents/Plink_1/plink

grep -v "^#" /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095.vcf | wc -l


#1. Descargas plink y en la misma carpeta pones tu vcf y corres la siguiente linea:

./plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 5 0.20 --out allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_prune_50_5_02

#2. PODAR EL ARCHIVO DE LOS SNPS CON LD:

./plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --extract allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_prune_50_5_02.prune.in --make-bed --pca --out /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_LD_prunned_50_5_02

#3.TRANSFORMAR EL PRUNEDDATA ( ARCHIVO PODADO) A VCF DE NUEVO

./plink2  --bfile /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_LD_prunned_50_5_02 --double-id --allow-extra-chr --set-missing-var-ids @:# --cow --recode vcf --out /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_LD_prunned_50_5_02

grep -v "^#" /Users/katherin_otalora/Documents/Plink_1/plink/allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS095_LD_prunned_50_5_02.vcf | wc -l


#Nos quedamos finalmente con el output de 50 5 0.20 con un total de 1785 SNPs para allSamples maxMISS0.95

#transform to .gz

bcftools view -O z -o allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.vcf.gz allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.vcf | wc -l


#..................................................................................................................................................................................................................................................................
#..................................................................................................................................................................................................................................................................

#check vcf sample names. Plink has doubled these ids

bcftools query -l /Users/Katherin_Otalora/Desktop/Katherin_PhD/Results/PLINK/plink/intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf

#or

bcftools view --header-only intSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS075.vcf | wc -l

#make a list of new single headers called new_names and rename

bcftools reheader -s New_names.txt -o new.bcf intSample.g.geno_filtered_LD.vcf
