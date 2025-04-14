# Population genetics and phylogeographic history of the insular lizard *Podarcis lilfordi* (Gunther, 1874) from the Balearic Islands based on genome-wide polymorphic data

This repository includes several codes for the identification of SNPs from GBS data, as well as the incorporation of SNPs from RAD-seq data reported by Bassitta et al. (2021). It also contains tools for downstream analyses, including genetic structure, phylogeographic patterns, diversity indices, and outlier detection. These codes were used to perform the analyses presented in the publication: https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.11407.



## SNP Discovery Workflow (GATK pipeline)


```mermaid
flowchart TD
    A["1. Reference genome assembly (Catalan Biogenome Project)"]
    A --> A1["18 autosomal chromosomes + sex chromosomes"]
    A1 --> A2["Female (ZW): Scaffold 14 = Z, Scaffold 20 = W"]
    A2 --> A3["Genome size: 1.4 Gb (Diploid)"]

    A3 --> B["2. Alignment of reads to genome"]
    B --> B1["Software: BWA, Samtools, Java, GATK"]

    B1 --> C["3. Variant Calling"]
    C --> C1["GATK HaplotypeCaller (--minimum-mapping-quality 20)"]
    C1 --> C2["CombineGVCFs"]
    C2 --> C3["GenotypeGVCFs"]
    C3 --> C4["Choose parameters depending on data"]

    C4 --> D["4. SNP Filtering"]
    D --> D1["VCFTOOLS v0.1.13, PLINK v1.9"]
    D1 --> D2["Filters: min-alleles, max-alleles, max-missing, minDP, maxDP, minQ, maf, indep-pairwise"]

    D2 --> E["5. Diversity Indexes"]
    E --> E1["STACKS v2.61, Hierfstat, PerformanceAnalytics, geodist (R 4.1.1)"]
    E1 --> E2["Indexes: Ar, Ho, He, Fis"]
    E2 --> E3["Pearson correlations with island variables"]
    E3 --> E4["Mantel test"]

    D2 --> F["6. Genetic Structuring"]
    F --> F1["Adegenet, dartR, ape (R 4.1.1)"]
    F1 --> F2["Methods: FST, PCA, K-means, DAPC"]
    F2 --> F3["Cross-validation: xvalDapc"]

    D2 --> G["7. Admixture Analysis"]
    G --> G1["ADMIXTURE v1.3.0"]
    G1 --> G2["Input: .bed file from PLINK"]

    D2 --> H["8. Contemporary Migration"]
    H --> H1["DiveMigrate"]
    H1 --> H2["Genepop format from PGDSpider"]

    D2 --> I["9. Phylogenetic Analysis"]
    I --> I1["Vcf2phylip.py, trimAI v1.4"]
    I1 --> I2["Geneious v2023.0.3, RAxML 8.2.12 (CIPRES), JModelTest 2.1.10"]
    I2 --> I3["IQ-TREE v2.2.0.8, Figtree v1.4.4"]

    D2 --> J["10. Detection of Outlier SNPs"]
    J --> J1["BayeScan v2.1 (1:100), PCAdapt (R)"]
    J1 --> J2["FST & log10(PO) filtering in R and Python"]
    J2 --> J3["Input: .gen from PGDSpider, then .txt"]

    J3 --> K["11. Genome-Environment Associations"]
    K --> K1["RDA (vegan package, R 4.1.1)"]

    J3 --> L["12. Enrichment Analysis of Outliers"]
    L --> L1["CDS annotation in R"]
    L1 --> L2["gGOSt in g:Profiler"]
    L2 --> L3["Visualization with REViGO"]
    L3 --> L4["GFF file from P. lilfordi genome"]

    D2 --> M["Datasets"]
    M --> M1["Combined: 6,394,354 SNPs (191 individuals)"]
    M1 --> M2["ExtSamples: 4,851,070 SNPs (91 individuals)"]
    M2 --> M3["IntSamples: 1,888,392 SNPs (100 individuals)"]
