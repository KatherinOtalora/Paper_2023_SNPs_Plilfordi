################### Identifying candidate genes #######################
# see script at https://github.com/rnnh/bioinfo-notebook/blob/master/scripts/annotating_snps.R

setwd("/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi")

##LOAD LIBRARIES 

library(dplyr)

##Define files and parameters

GFF_file <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/PODLIA.FA.gff3"
SNP_file <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/names_outliers_filtred_LOGPO_FST_extamples_withoutLD_1_100_K=10.txt"
#upstream.int <- as.integer("25000") #The variable 'upstream.int' is used to determine how far upstream from an
#     annotated feature a SNP can be. This can be set to 0 if you do not want
#     upstream SNPs to be considered. Setting it to 1000 will mean that SNPs
#     up to 1,000 bases/1kb upstream from a feature will be annotated.
#downstream.int <- as.integer("1000")


# Reading assembly annotation file (.gff)
gff <- read.delim(
  GFF_file,
  header = T, sep="\t"
)
colnames(gff) <- c("chr", "source", "feature", "start", "end", "score",
                   "strand", "frame", "attribute")
#Miramos cuandos tipos de genes CDS exon o transcript hay en mi gff

table(gff$feature)

# Reading SNPs 
SNPs.df <- read.delim(
  SNP_file,
  header = T, sep="\t"
)
colnames(SNPs.df) <- c("chr", "POS")


# Joining data frames using dplyr ==============================================

# In this section...
#   - inner_join() is used to join together the genome annotation and SNPs
#     data frames along the column "chr": i.e. rows with the same value
#     for "chr" are joined together
#   - filter() is used to filter out SNPs which do not fall within regions of
#     interest: i.e. SNPs that are not within- or upstream of- features

# Joining data frames with genome annotation and SNPs
SNPs_with_annotations.df <- inner_join(gff, SNPs.df, 
                                       by = "chr") %>%
  filter(POS >=start &
           POS <= end)
#filter(POS >= (start - upstream.int) &
#       POS <= (end + downstream.int))

#Te saca un fichero de output donde en la última columna (POS) hay la posición de lo SNP (outlier o simplemente cualquier SNP te interese ver) 
#Y la anotación de lo que está cerca adentro de la ventana (hay varias filas por cada SNP/POS).
#Hay varias filas por POS o SNP porque se puede encontrar en exones CDS y genes, por eso ahora hay que filtrar por genes o lo que nos interese


# Ordering filtered data frame of SNPs with annotations ========================
#La función attach() en R Language se utiliza para acceder a las variables presentes en el marco de datos sin llamar al marco de datos
attach(SNPs_with_annotations.df)
SNPs_with_annotations.df <- SNPs_with_annotations.df[order(chr, start, end), ]
detach(SNPs_with_annotations.df)


# Exporting SNPs with general annotations to tab-separated value (.tsv) file ===========
write.table(SNPs_with_annotations.df,
            file = "Output_extsamples_anotacion.tsv",
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)

#select CDS only , CDS = gen - intrones - UTRs , es decir solo regiones que codifican a proteina

SNPs_with_annotations_genes.df  <- SNPs_with_annotations.df  %>% filter(feature == "gene")
SNPs_with_annotations_genes.df 

SNPs_with_annotations_CDS.df  <- SNPs_with_annotations.df  %>% filter(feature == "CDS")
SNPs_with_annotations_CDS.df 

output_name_1 <- "Output_extsamples_anotacion_filtredgenes.tsv" 
output_name_2 <- "Output_extsamples_anotacion_filtredCDS.tsv"

# Exporting SNPs with annotations GENES to tab-separated value (.tsv) file ===========
write.table(SNPs_with_annotations_genes.df,
            file = output_name_1,
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)

# Exporting SNPs with annotations CDS to tab-separated value (.tsv) file ===========
write.table(SNPs_with_annotations_CDS.df,
            file = output_name_2,
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)

#Al outut que salio la columna attribute le elimino las palabras ID= y ;gene_type=protein_coding para quedarme solo con el nombre de la proteina 

# Removing redundant data frames
rm(gff, SNPs.df)


....................................................................................................

###################################### DEINIFICION DE PROTEINAS CDS #####################################
#Vamos a hacer merge de los SNPs outliers genes con el PODLIA.protein_definitions para sacar la informacion de codificacion de la proteina, Pegar las  columnas con las informacion
#Iniciamos leyendo nuestro archivo editado donde eliminamos el ID= y quedamos solo con el nombre de la proteina

SNP_file_1 <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/Output_extsamples_anotacion_filtredCDS.tsv"

# Reading SNPs con su anotacion del gff y editados (sin ID= y ;gene_type=protein_coding)
SNPs_1.df <- read.delim(
  SNP_file_1,
  header = T, sep="\t"
)
#No es necesario nombrar los colnames

# Reading Definicion de las proteinas ( se le quito el P1 o P2,3,4,5,6,7,8,9 del final del nombre de cada proteina para el merge en .txt)
Protein_def_file <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/PODLIA.protein_definitions_2.txt"
Protein_def_file <- read.delim(
  Protein_def_file,
  header = T, sep="\t"
)
colnames(Protein_def_file)=c("attribute","protein_def")

#Usamos merge para pegar las dos variables 
x=SNPs_1.df
y=Protein_def_file 
z=merge(x,y, by="attribute")
head (z)

# Exporting SNPs with annotations to tab-separated value (.tsv) file ===========
write.table(z,
            file = "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/CDS_annotation_outliers_filtred_LOGPO_FST_extsamples_withoutLD_1_100_K=10.tsv",
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)
#Al outut que salio la columna attribute le elimino las palabras ID= y ;gene_type=protein_coding para quedarme solo con el nombre de la proteina 

....................................................................................................

###################################### DEINIFICION DE PROTEINAS GENE #####################################
#Vamos a hacer merge de los SNPs outliers genes con el PODLIA.protein_definitions para sacar la informacion de codificacion de la proteina, Pegar las  columnas con las informacion
#Iniciamos leyendo nuestro archivo editado donde eliminamos el ID= y quedamos solo con el nombre de la proteina

SNP_file_1 <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/Output_extsamples_anotacion_filtredgenes.tsv"

# Reading SNPs con su anotacion del gff y editados (sin ID= y ;gene_type=protein_coding)
SNPs_1.df <- read.delim(
  SNP_file_1,
  header = T, sep="\t"
)
#No es necesario nombrar los colnames

# Reading Definicion de las proteinas ( se le quito el P1 o P2,3,4,5,6,7,8,9 del final del nombre de cada proteina para el merge en .txt)
Protein_def_file <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/PODLIA.protein_definitions_2.txt"
Protein_def_file <- read.delim(
  Protein_def_file,
  header = T, sep="\t"
)
colnames(Protein_def_file)=c("attribute","protein_def")

#Usamos merge para pegar las dos variables 
x=SNPs_1.df
y=Protein_def_file 
z=merge(x,y, by="attribute")
head (z)

# Exporting SNPs with annotations to tab-separated value (.tsv) file ===========
write.table(z,
            file = "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/GENE_annotation_outliers_filtred_LOGPO_FST_extsamples_withoutLD_1_100_K=10.tsv",
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)
#Al outut que salio la columna attribute le elimino las palabras ID= y ;gene_type=protein_coding para quedarme solo con el nombre de la proteina 

....................................................................................................
################################ Anotacion con LncRNA gff #########################
# Reading Definicion de  los Long_non_coding_RNAs( se le quito el ### ENTRE DATOS para  que se pueda hacer LA ANOTACION)
Long_non_coding_RNAs_gff_file <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/PODLIncA_2.gff3"
LncRNA_gff_file <- read.delim(
  Long_non_coding_RNAs_gff_file,
  header = T, sep="\t"
)
colnames(LncRNA_gff_file ) <- c("chr", "source", "feature", "start", "end", "score",
                                "strand", "frame", "attribute")

#upstream.int.lncRNA <- as.integer("25000") #The variable 'upstream.int' is used to determine how far upstream from an
#     annotated feature a SNP can be. This can be set to 0 if you do not want
#     upstream SNPs to be considered. Setting it to 1000 will mean that SNPs
#     up to 1,000 bases/1kb upstream from a feature will be annotated.

# Joining data frames with genome annotation and SNPs
SNPs_annotations_LncRNA.df <- inner_join(LncRNA_gff_file, SNPs.df, 
                                         by = "chr") %>%
  filter(POS >= start &
           POS <= end)
#  filter(POS >= (start - upstream.int.lncRNA) &
#           POS <= end)
# Ordering filtered data frame of SNPs with annotations ========================
#La función attach() en R Language se utiliza para acceder a las variables presentes en el marco de datos sin llamar al marco de datos
attach(SNPs_annotations_LncRNA.df)
SNPs_annotations_LncRNA.df <- SNPs_annotations_LncRNA.df[order(chr, start, end), ]
detach(SNPs_annotations_LncRNA.df)

#select genes only
#SNPs_annotations_genes_LncRNA.df  <- SNPs_annotations_LncRNA.df %>% filter(feature == "gene")
#SNPs_annotations_genes_LncRNA.df 

output_name_3 <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/Output_extsamples_anotacion_lncRNA.tsv" 

# Exporting SNPs with annotations to tab-separated value (.tsv) file ===========
write.table(SNPs_annotations_LncRNA.df ,
            file = output_name_3,
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)
....................................................................................................

######################################ANOTACION DE GO CDS #####################################
#Vamos a hacer merge de las Proteinas para tener el GO anotacion
#Iniciamos leyendo nuestro archivo editado donde eliminamos y usaremos la columna que tiene el nombre de la proteina

SNP_file_protein <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/CDS_annotation_outliers_filtred_LOGPO_FST_extsamples_withoutLD_1_100_K=10.tsv"

# Reading SNPs con su anotacion del gff y editados (sin ID= y ;gene_type=protein_coding)
SNP_file_protein.df <- read.delim(
  SNP_file_protein,
  header = T, sep="\t"
)
#No es necesario nombrar los colnames

# Reading Definicion deL GO file ( se le quito el P1 o P2,3,4,5,6,7,8,9 del final del nombre de cada proteina para el merge en .txt)
GO_def_file <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/PODLIA.GO_prot_1.txt"
GO_def_file_1 <- read.delim(
  GO_def_file,
  header = T, sep="\t"
)
colnames(GO_def_file_1)=c("attribute","GO_def")

#Usamos merge para pegar las dos variables 
x=SNP_file_protein.df
y=GO_def_file_1
z=merge(x,y, by="attribute")
head (z)

# Exporting SNPs with annotations to tab-separated value (.tsv) file ===========
write.table(z,
            file = "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/CDS_annotation_GO_outliers_filtred_LOGPO_FST_extsamples_withoutLD_1_100_K=10.tsv",
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)


######################################ANOTACION DE GO GENE #####################################
#Vamos a hacer merge de las Proteinas para tener el GO anotacion
#Iniciamos leyendo nuestro archivo editado donde eliminamos y usaremos la columna que tiene el nombre de la proteina

SNP_file_protein <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/GENE_annotation_outliers_filtred_LOGPO_FST_extsamples_withoutLD_1_100_K=10.tsv"

# Reading SNPs con su anotacion del gff y editados (sin ID= y ;gene_type=protein_coding)
SNP_file_protein.df <- read.delim(
  SNP_file_protein,
  header = T, sep="\t"
)
#No es necesario nombrar los colnames

# Reading Definicion deL GO file ( se le quito el P1 o P2,3,4,5,6,7,8,9 del final del nombre de cada proteina para el merge en .txt)
GO_def_file <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/PODLIA.GO_prot_1.txt"
GO_def_file_1 <- read.delim(
  GO_def_file,
  header = T, sep="\t"
)
colnames(GO_def_file_1)=c("attribute","GO_def")

#Usamos merge para pegar las dos variables 
x=SNP_file_protein.df
y=GO_def_file_1
z=merge(x,y, by="attribute")
head (z)

# Exporting SNPs with annotations to tab-separated value (.tsv) file ===========
write.table(z,
            file = "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/GENE_annotation_GO_outliers_filtred_LOGPO_FST_extsamples_withoutLD_1_100_K=10.tsv",
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)


##################### Antes de hacer este paso hay que hacer el analisis de enriquecimiento en la pagina #################
#### Se debe sacar el output con los term id y los intersectos mas significativos que seran nuestro GO_def_file ###################
###################################### DEINIFICION DE GO CDS #####################################

SNP_file_GO <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/CDS_annotation_GO_outliers_filtred_LOGPO_FST_extsamples_withoutLD_1_100_K=10.tsv"

# Reading SNPs con su anotacion del gff y editados (sin ID= y ;gene_type=protein_coding)
SNPs_1.df <- read.delim(
  SNP_file_GO,
  header = T, sep="\t"
)
#No es necesario nombrar los colnames

# Reading Definicion de las proteinas ( se le quito el P1 o P2,3,4,5,6,7,8,9 del final del nombre de cada proteina para el merge en .txt)
GO_def_file <- "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/gProfiler_pmuralis_06-07-2023_12-10-33__intersections_extsamples_onpoint.csv"
GO_def_file <- read.delim(
  GO_def_file,
  header = T, sep=","
)
colnames(SNPs_1.df) <- c("attribute","chr", "source", "feature", "start", "end", "score",
                         "strand", "frame", "POS", "protein_def","term_id")

#Usamos merge para pegar las dos variables 
x=SNPs_1.df
y=GO_def_file
z=merge(x,y, by="term_id")
head (z)

# Exporting SNPs with annotations to tab-separated value (.tsv) file ===========
write.table(z,
            file = "/Users/katherin_otalora/Documents/GFF_Genome_Plilfordi/Anotacion_extsamples/CDS_annotation_term_id_GO_outliers_filtred_LOGPO_FST_extsamples_withoutLD_1_100_K=10.tsv",
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)


# Exiting ======================================================================
quit(save = "no")
