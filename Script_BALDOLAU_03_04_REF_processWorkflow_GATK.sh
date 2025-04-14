WDIR=$(pwd)  #Declarar la variable WDIR con el pwd ( es decir donde se encuentra ubicado), para mas adelante ponerlo como directorio y guardar los datos en un control flow

# Get FASTA #Entrar al genoma de P.lilfordi .fa
#Para crear enlaces utilizamos el comando ln y la opción -s para especificar enlaces simbólicos.
ln -s /scratch/devel/talioto/denovo_assemblies/podarcis_lilfordi/assemblies/rPodLil1.1/rPodLil1.1.scaffolds.fa .

#Se intenta vincular el directorio donde esta el genoma con a la carpeta actual (por eso se pone un punto al final), punto signfica hacerlo aqui.

# Create FASTA index and dict
#SamTools es un conjunto de utilidades para interactuar y postprocesar alineamientos de lectura de secuencias de ADN cortas en los formatos SAM, BAM y CRAM
samtools faidx rPodLil1.1.scaffolds.fa #indexa o consulta regiones de un archivo fasta
samtools dict rPodLil1.1.scaffolds.fa > rPodLil1.1.scaffolds.dict #crea un archivo o output de diccionario de secuencias a partir de un archivo fasta donde te da los datos de etiquetas, de las secuencias, incluso de los header

#Dict tiene varias opciones que son:

#Especificar el nombre del conjunto para la etiqueta AS.
#Especificar el nombre de la especie para la etiqueta SP.
#Utilizar el alias o nombre alternativo para las secuencias.
#Imprimir o no la línea de header @HD.

# Create BWA index
bwa index rPodLil1.1.scaffolds.fa #Todos los alineadores NGS necesitan que las secuencias de referencia estén indexadas. Usted utilizaría el nombre base del índice con sus reads fastq en el momento del aleamiento

# Get external FASTQs https://www.ebi.ac.uk/ena/browser/view/PRJNA645796 > ena_files.zip
unzip ena_files.zip #Descromprimir el archivo ena_files.zip que son los datos de Bassita etal bajo el proyecto PRJNA645796
mkdir extSample #Hacer una carpeta de extSample
mv SRR122* extSample/. #Mover todos los archivos que tengan SRR122 a la carpeta extSample

# Create sample sheet (one line per sample ID)
intSampleSheet.tsv #Crear una hoja con valores separados por tabuladores para el intSample and extSample
extSampleSheet.tsv

# Get internal FASTQs
#Se rellenan el archivo intSample.tsv con la informacion de los barcode que aparecen en los nombres de cada fastq.gz de los data set

while read barcode; do \
mkdir -p intSample/$barcode; \
cat /scratch/devel/jrtrotta/BALDOLAU/BALDOLAU_03_04/BALDOLAU_03_04/*/*_fastqs/${barcode}_R1_.fastq.gz > intSample/$barcode/${barcode}_1.fastq.gz; \
cat /scratch/devel/jrtrotta/BALDOLAU/BALDOLAU_03_04/BALDOLAU_03_04/*/*_fastqs/${barcode}_R2_.fastq.gz > intSample/$barcode/${barcode}_2.fastq.gz; \
done < intSampleSheet.tsv


# Align FASTQS
#Alineamiento de los fastq con el genoma de referencia
#Primero define sheet que seria cada archivo .tsv mientras lee el barcode entra a la carpeta barcode de cada tipo de muestra
#Segundo usa BWA-MEM que es adecuado para alinear lecturas largas de alta calidad que van desde 70 bp a 1 Mbps contra un gran genoma de referencia
# Luego se crea una tuberia con | samtools para postprocesar los alineamientos  creando un archivo .bam y luego usando Java
#Luego se usa samtools index, para indexar el archivo BAM ordenado y asi extraer rápidamente los alineamientos que se solapan con determinadas regiones genómicas. Además, la indexación es necesaria para los visores del genoma, como IGV, para que los visores puedan mostrar rápidamente los alineamientos en cada región genómica hacia la que se navega.
#Luego remueve los .tmp.bam
#Luego se ubica en la carpeta WDIR inicialmente declarada
#Luego done para parar el loop

for sheet in $(ls *SampleSheet.tsv); do \
while read barcode; do \
cd ${sheet%Sheet.tsv}/$barcode; \
bwa mem -t 4 -M -R "@RG\tID:$barcode\tSM:$barcode\tLB:$barcode\tPL:ILLUMINA\tPU:HISEQ" \
/scratch/devel/jrtrotta/BALDOLAU/BALDOLAU_03_04_REF/rPodLil1.1.scaffolds.fasta \
${barcode}_1.fastq.gz ${barcode}_2.fastq.gz | samtools view -S -b -> $barcode.tmp.bam; \
java -Xmx60g -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/1.110/SortSam.jar I=$barcode.tmp.bam O=$barcode.bam SO=coordinate; \
samtools index $barcode.bam; \
rm $barcode.tmp.bam; \
cd $WDIR; \
done < $sheet; \
done

# Haplotype Caller
#El HaplotypeCaller es capaz de llamar SNPs e indels simultáneamente a través del ensamblaje local de-novo de haplotipos en una región activa.
#En otras palabras, cada vez que el programa encuentra una región que muestra signos de variación, descarta la información de mapeo existente y vuelve a ensamblar completamente las lecturas en esa región.
#Esto permite al HaplotypeCaller ser más preciso cuando llama a regiones que tradicionalmente son difíciles de llamar, por ejemplo, cuando contienen diferentes tipos de variantes cerca unas de otras. También hace que el HaplotypeCaller sea mucho mejor a la hora de determinar las indels

for sheet in $(ls *SampleSheet.tsv); do \
while read barcode; do \
cd ${sheet%Sheet.tsv}/$barcode; \
gatk --java-options "-Xmx38g -Xms38g" HaplotypeCaller \
-input $barcode.bam \
--reference /scratch/devel/jrtrotta/BALDOLAU/BALDOLAU_03_04_REF/rPodLil1.1.scaffolds.fasta \
--emit-ref-confidence GVCF \
--minimum-mapping-quality 20 \
--output $barcode.g.vcf.gz
cd $WDIR; \
done < $sheet; \
done

# Combine gVCF (BALDOLAU+extData 1 BALDOLAU only)
# gVCF son las siglas de Genomic VCF. Un GVCF es un tipo de VCF, por lo que la especificación básica del formato es la misma que la de un VCF normal, pero un VCF Genómico contiene información extra.
# La diferencia clave entre un VCF normal y un GVCF es que el GVCF tiene registros para todos los sitios, haya o no una llamada de variante allí. El objetivo es tener todos los sitios representados en el archivo para poder hacer un análisis conjunto de una cohorte en los pasos siguientes. Los registros de un GVCF incluyen una estimación precisa de la confianza que tenemos en la determinación de que los sitios son homocigotos de referencia o no.
# Combine gVCF : Está pensado para ser utilizado para la fusión de GVCFs que eventualmente serán introducidos en GenotypeGVCFs.
# Se combinan el genoma de referencia con todos los outputs g.vcf.gz del haplotype caller (allSamples)
# Se combinan el genoma de referencia con todos los outputs g.vcf.gz del haplotype caller (intSamples)

allVariants=$(for sample in $(ls /scratch/devel/jrtrotta/BALDOLAU/BALDOLAU_03_04_REF/*Sample/*/*.g.vcf.gz); do echo $sample; done | awk '{print "--variant "$0" "}')
gatk --java-options "-Xmx56g" CombineGVCFs \
--reference /scratch/devel/jrtrotta/BALDOLAU/BALDOLAU_03_04_REF/rPodLil1.1.scaffolds.fasta \
--output allSample.g.vcf.gz \
$allVariants

intVariants=$(for sample in $(ls /scratch/devel/jrtrotta/BALDOLAU/BALDOLAU_03_04_REF/intSample/*/*.g.vcf.gz); do echo $sample; done | awk '{print "--variant "$0" "}')
gatk --java-options "-Xmx56g" CombineGVCFs \
--reference /scratch/devel/jrtrotta/BALDOLAU/BALDOLAU_03_04_REF/rPodLil1.1.scaffolds.fasta \
--output intSample.g.vcf.gz \
$intVariants

# Genotype gVCF (BALDOLAU+extData 1 BALDOLAU only)
#Realiza un genotipado en una o más muestras prellamadas con HaplotypeCaller,
# Genotype gVCF: utiliza las variantes potenciales del HaplotypeCaller y realiza el genotipado conjunto.
#Mirará la información disponible para cada sitio de los alelos variantes y no variantes a través de todas las muestras, y producirá un archivo VCF que contiene sólo los sitios que encontró como variantes en al menos una muestra.
#Utiliza el genoma de referencia, el nombre del output que tendra unicamente los sitios que encontro con variantes en al menos una muestras
#Utiliza el --variant que es el output del Haplotypecaller que son las variantes potenciales

gatk --java-options "-Xmx56g" GenotypeGVCFs \
--reference /scratch/devel/jrtrotta/BALDOLAU/BALDOLAU_03_04_REF/rPodLil1.1.scaffolds.fasta \
--output allSample.g.geno.vcf.gz \
--variant allSample.g.vcf.gz

gatk --java-options "-Xmx56g" GenotypeGVCFs \
--reference /scratch/devel/jrtrotta/BALDOLAU/BALDOLAU_03_04_REF/rPodLil1.1.scaffolds.fasta \
--output intSample.g.geno.vcf.gz \
--variant intSample.g.vcf.gz

# Get supported variants
#Se usa esta herramienta gerSupported para hacer un filtrado de los vcf con -D que es la profundidad o deep de la secuenciacion del SNPs Este es esencialmente el número de reads que se han asignado a esta posición.
#-N Numero de muestras que pasan los filtros y se aceptan las variantes
#-I Es el input en este caso el vcf al que quieres hacer el filtrado
#-O es el nombre del ouput donde quieres que se guarde la informacion filtrada
/scratch/production/DAT/apps/GETSUPPORTEDVARIANTS/getSupportedVariants -D 10 -N 20 -I intSample.g.geno.vcf.gz -O intSample.g.geno.supp.DP10.altDP2.AF005.n20.vcf.gz
/scratch/production/DAT/apps/GETSUPPORTEDVARIANTS/getSupportedVariants -D 10 -N 38 -I allSample.g.geno.vcf.gz -O allSample.g.geno.supp.DP10.altDP2.AF005.n38.vcf.gz

#usage: getSupportedVariants [options]

# Utility to select variants in which at least S samples with a genotype not
# equal to ./. or 0/0 are supported by a depth of D, GQ>=G with A reads different
# to the reference allele and a frequency >=F

#optional arguments:
#  -h, --help            show this help message and exit
#  --minAltDepth MINALTDEPTH, -A MINALTDEPTH
#                        INT Minimum sample alternative allele depth [2]
#  --minAltFrequency MINALTFREQUENCY, -F MINALTFREQUENCY
#                        INT Minimum sample alternative allele frequency [0.05]
#  --minDepth MINDEPTH, -D MINDEPTH
#                        INT Minimum sample read depth [10]
#  --minGenotypeQuality MINGENOTYPEQUALITY, -G MINGENOTYPEQUALITY
#                        INT Minimum GQ [0]
#  --numberSamples NUMBERSAMPLES, -N NUMBERSAMPLES
#                        INT Number of samples passing filters to accept the
#                        variant [1]
#  --tag, -t             BOOL Tag variants not passing this filter instead of
#                        removing them
#  --vcfDiscarded VCFDISCARDED, -d VCFDISCARDED
#                        FILE Path to the vcf file to save discarded variants
#                        (debug option)
#  --vcfIn VCFIN, -I VCFIN
#                        FILE Path to the vcf file to classify
#  --vcfOut VCFOUT, -O VCFOUT
#                        FILE Name of the output vcf file
