######### IQTREE ############

brew install boost
brew install eigen
brew install libomp


git clone https://github.com/iqtree/iqtree2.git

git clone --recursive https://github.com/Cibiv/IQ-TREE.git
cd IQ-TREE
git checkout latest
git submodule init
git submodule update

/Applications/CMake.app/Contents/bin/cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ .

# Se descarga de http://www.iqtree.org/#download

cd /Users/katherin_otalora/Documents/iqtree-2.2.0.8-MacOSX/bin/iqtree2
/Users/katherin_otalora/Documents/iqtree-2.2.0.8-MacOSX/bin/iqtree2 --help

################# Inferring rooted trees without outgroup #####################

################### PRUEBA CON DATOS DEL MANUAL #####################

#Descargamos los archivos y corrimos el analisis

cd /Users/katherin_otalora/Documents/iqtree-2.2.0.8-MacOSX/bin

./iqtree2 -s bovidae.phy -p bovidae.nex -B 1000 -T AUTO -pre rev_dna_bovidae

./iqtree2 -s bovidae.phy -p rev_dna_bovidae.best_scheme.nex --model-joint 12.12 -B 1000 -T AUTO -pre nonrev_dna_prueba

#Efectivamente salio el output con el .rootstrap

########## ANALISIS CON NUESTRO DATA SET CONSENSUS ###########################

#Inferir arbol no enraizado con modelos NO reversibles, non-reversible DNA model

./iqtree2  -s allsamples_consensus_15samples_gt096_plilfordi.fasta --model-joint 12.12 -B 1000 -T AUTO -pre NOrev_allsamples_B1K

# Probaremos de nuevo inferir el arbol pero con un Bootstrap de 10000 y 100000

./iqtree2  -s allsamples_consensus_15samples_gt096_plilfordi.fasta --model-joint 12.12 -B 10000 -T AUTO -pre NOrev_allsamples_B10K
./iqtree2  -s allsamples_consensus_15samples_gt096_plilfordi.fasta --model-joint 12.12 -B 100000 -T AUTO -pre NOrev_allsamples_B100K

# Vamos a intentar el análisis eliminando los sitios invariantes con el siguiente script, recuerde que solo puedo sacar el variantes con el IQTREE. 2.2.0.8:

./iqtree2  -s allsamples_consensus_15samples_gt096_plilfordi.fasta -m GTR+ASC -B 1000 -T AUTO -pre allsamples_consensus_15samples_gt096_plilfordi
#ERROR y me saca un allsamples_consensus_15samples_gt096_plilfordi.varasites.phy con los sitios invariantes

#Inferir arbol no enraizado con modelos NO reversibles, non-reversible DNA model

./iqtree2  -s allsamples_consensus_15samples_gt096_plilfordi.varasites.phy --model-joint 12.12 -B 1000 -T AUTO -pre NOrev_allsamples_B1K

# Probaremos de nuevo inferir el arbol pero con un Bootstrap de 10000 y 100000

./iqtree2  -s allsamples_consensus_15samples_gt096_plilfordi.varasites.phy --model-joint 12.12 -B 10000 -T AUTO -pre NOrev_allsamples_B10K
./iqtree2  -s allsamples_consensus_15samples_gt096_plilfordi.varasites.phy --model-joint 12.12 -B 100000 -T AUTO -pre NOrev_allsamples_B100K

########## ANALISIS CON NUESTRO DATA SET COMPLETO DE TODOS LOS SNPS ###########################

#Inferir arbol no enraizado con modelos NO reversibles, non-reversible DNA model

./iqtree2  -s 1_allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_wtoutsamplescomplot_curated_183samples_gt096_N.min4.fasta --model-joint 12.12 -B 1000 -T AUTO -pre NOrev_1_allsamples_complete_B1K

# Error al correr el análisis por esto pondré UNREST y no 12.12

./iqtree2  -s 1_allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_wtoutsamplescomplot_curated_183samples_gt096_N.min4.fasta --model-joint UNREST -B 1000 -T AUTO -pre NOrev_1_allsamples_complete_B1K

# Probaremos de nuevo inferir el arbol pero con un Bootstrap de 10000 y 100000

./iqtree2  -s 1_allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_wtoutsamplescomplot_curated_183samples_gt096_N.min4.fasta --model-joint UNREST -B 10000 -T AUTO -pre NOrev_1_allsamples_complete_B10K

./iqtree2  -s 1_allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_wtoutsamplescomplot_curated_183samples_gt096_N.min4.fasta --model-joint UNREST -B 100000 -T AUTO -pre NOrev_1_allsamples_complete_B100K

# Vamos a intentar el análisis eliminando los sitios invariantes con el siguiente script, recuerde que solo puedo sacar el variantes con el IQTREE. 2.2.0.8:

./iqtree2  -s 1_allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_wtoutsamplescomplot_curated_183samples_gt096_N.min4.varsites.phy --model-joint UNREST -B 1000 -T AUTO -pre NOrev_1_allsamples_varsites_B1K

# Probaremos de nuevo inferir el arbol pero con un Bootstrap de 10000 y 100000

./iqtree2  -s 1_allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_wtoutsamplescomplot_curated_183samples_gt096_N.min4.varsites.phy --model-joint UNREST -B 10000 -T AUTO -pre NOrev_1_allsamples_varsites_B10K

./iqtree2  -s 1_allSample.plilfordi_filter_minalle2_maxalle2_MAF005_Q30_minDP10_maxD50_maxMISS090_LD_prunned_50_5_02_wtoutsamplescomplot_curated_183samples_gt096_N.min4.varsites.phy --model-joint UNREST -B 100000 -T AUTO -pre NOrev_1_allsamples_varsites_B100K
