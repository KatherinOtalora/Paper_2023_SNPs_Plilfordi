###First of all, remove individuals with more than xx % of Ns.
#Vamos a la carpeta donde tengo mi informacion

cd  /Users/katherin_otalora/Documents/Tree_SNPs/Tree_SNPs_Script

#First Calculate proportion of missing data per individual

sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.vcf > allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_v3.vcf
vcftools --gzvcf allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_v3.vcf.gz --missing-indv --out allsamples_missing_ind #it generates a .imiss output

##load the .imiss on excel and check the Fmiss (fraction of missing data). There are 8 individuals with more than 10% Ns
AR6515
AR6519
AR6482
AR6517
AR6516
SRR12207018
SRR12207035
SRR12207056

#Remove 8 specimens with more than 10% Ns
bcftools view -s ^AR6515_AR6515,AR6519_AR6519,AR6482_AR6482,AR6517_AR6517,AR6516_AR6516,SRR12207018_SRR12207018,SRR12207035_SRR12207035,SRR12207056_SRR12207056 allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02.vcf.gz > allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8filtered.vcf #use ^ to exclude samples
bcftools query -l allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8filtered.vcf

#extract the outliers from allsamples
cd  /Users/katherin_otalora/Documents/Tree_SNPs/Tree_SNPs_Script
#PARA Pasar a la version actual v4.3
sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8filtered.vcf > allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8filtered_V3.vcf
grep  -c -v "^#" allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8filtered_V3.vcf
#2876
#vcftools --gzvcf allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8filtered_V3.vcf --exclude-positions 88655044 > allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8_outlier_filtered.vcf
bcftools view -T ^ls_exclude_outliers_allsamples_LD.txt allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8filtered_V3.vcf > allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8_outlier_filtered.vcf
grep  -c -v "^#" allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8_outlier_filtered.vcf
#2867 ELIMINO SOLO 9 de 18 PORQUE AL HACER EL LD ya habia eliminado los otros 9 outliers descritos para allsamples

##see https://gitee.com/khjia/vcf2phylip
# transform vcf to phylyp. The probelm here it is that this transformation does not consider indels, it just replace them with Ns

cd  /Users/katherin_otalora/Documents/Tree_SNPs/Tree_SNPs_Script
python ./vcf2phylip.py  -i allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8_outlier_filtered.vcf -o allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8_outlier_filtered.phy

## Load the file in GENEIOUS to check the alignment. Save it to fasta

#Replace N with -

sed 's/N/-/g' allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8_outlier_filtered.min4.phy > allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8_outlier_filtered_gaps.phy


##using trimal, remove sites with gaps
cd /Users/katherin_otalora/Documents/Tree_SNPs/trimal-trimAl/source
#it removes columns with some percentage of gaps
./trimal -in allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_N8_outlier_filtered_gaps.phy -gt 0.95 -out all_nogaps095_filtredLD.fasta  -fasta
#2009 SNPS

#Replace - with N
sed 's/-/N/g' all_nogaps095_filtredLD.fasta > all_nogaps095_filtredLD_final.fasta
#2009 SNPS

#Las secuencias consenso se hicieron en Geneious


################El siguiente script no se utiló en el analisis #########
#Create a fasta for each population. Make a list of ids in txt file
cd /Users/katherin_otalora/Documents/Tree_SNPs/Tree_SNPs_Script/seqtk

./seqtk subseq all_nogaps095_filtredLD_final_names.fasta adgr.txt > adgr.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta aire.txt > aire.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta rei.txt > rei.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta dragonera.txt > dragonera.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta cabrera.txt > cabrera.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta colomer.txt > colomer.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta colom.txt > colom.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta porros.txt > porros.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta prcv.txt > prcv.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta ravl.txt > ravl.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta foradada.txt > foradada.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta escla.txt > escla.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta moltona.txt > moltona.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta guardia.txt > guardia.fas
./seqtk subseq all_nogaps095_filtredLD_final_names.fasta encurt.txt > encurt.fas


#### R ###
##Change sequence IDs, attaching islet code to fasta files

setwd("/Users/katherin_otalora/Documents/Tree_SNPs/Tree_SNPs_Script/seqtk")
install.packages("Biostrings")
library(Biostrings)

  myseq=readDNAStringSet("adgr.fas", "fasta") # im
  newids=gsub(" 2015", "_ADGR", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "adgr_aln.fas", format="fasta")

  myseq=readDNAStringSet("aire.fas", "fasta") # im
  newids=gsub(" 2015", "_aire", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "aire_aln.fas", format="fasta")

  myseq=readDNAStringSet("adgr.fas", "fasta") # im
  newids=gsub(" 2015", "_ADGR", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "adgr_aln.fas", format="fasta")

  myseq=readDNAStringSet("colom.fas", "fasta") # im
  newids=gsub(" 2015", "_colom", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "colom_aln.fas", format="fasta")

  myseq=readDNAStringSet("colomer.fas", "fasta") # im
  newids=gsub(" 2015", "_colomer", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "colomer_aln.fas", format="fasta")

  myseq=readDNAStringSet("foradada.fas", "fasta") # im
  newids=gsub(" 2015", "_foradada", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "foradada_aln.fas", format="fasta")

  myseq=readDNAStringSet("ravl.fas", "fasta") # im
  newids=gsub(" 2015", "_ravl", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "ravl_aln.fas", format="fasta")

  myseq=readDNAStringSet("prcv.fas", "fasta") # im
  newids=gsub(" 2015", "_prcv", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "prcv_aln.fas", format="fasta")

  myseq=readDNAStringSet("guardia.fas", "fasta") # im
  newids=gsub(" 2015", "_guardia", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "guardia_aln.fas", format="fasta")

  myseq=readDNAStringSet("moltona.fas", "fasta") # im
  newids=gsub(" 2015", "_moltona", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "moltona_aln.fas", format="fasta")


  myseq=readDNAStringSet("rei.fas", "fasta") # im
  newids=gsub(" 2015", "_rei", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "rei_aln.fas", format="fasta")

  myseq=readDNAStringSet("escla.fas", "fasta") # im
  newids=gsub(" 2015", "_escla", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "escla_aln.fas", format="fasta")


  myseq=readDNAStringSet("encurt.fas", "fasta") # im
  newids=gsub(" 2015", "_encurt", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "encurt_aln.fas", format="fasta")

  myseq=readDNAStringSet("dragonera.fas", "fasta") # im
  newids=gsub(" 2015", "_dragonera", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "dragonera_aln.fas", format="fasta")

  myseq=readDNAStringSet("cabrera.fas", "fasta") # im
  newids=gsub(" 2015", "_cabrera", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "cabrera_aln.fas", format="fasta")

  myseq=readDNAStringSet("porros.fas", "fasta") # im
  newids=gsub(" 2015", "_porros", names(myseq)) #to substitute a name with nothing (so to remove it). Use the number corresponding to the length of the aln
  names(myseq)=newids
  writeXStringSet(myseq, "porros_aln.fas", format="fasta")

  #Make final aln with all sequences with their islet code
  cat *_aln.fas> All_aln_95gaps.fas

#Make the consensus sequence for each population

setwd("~/Dropbox/RESEARCH/PODARCIS/Podarcis_SNPs/SNPs_data_analysis/Phylogeny/Islet_fasta")
library(Biostrings)
?consensusString

files=list.files(pattern=".fas$")
for (i in files) {
myseq=readDNAStringSet(i, "fasta") # im
cons=consensusString(myseq)
assign(i, cons)
write.table(cons, paste(i, c(".out"), sep=""), quote=FALSE, sep="\t")
}
