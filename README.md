# Population genetics and phylogeographic history of the insular lizard *Podarcis lilfordi* (Gunther, 1874) from the Balearic Islands based on genome-wide polymorphic data

This repository includes several codes for the identification of SNPs from GBS data, as well as the incorporation of SNPs from RAD-seq data reported by Bassitta et al. (2021). It also contains tools for downstream analyses, including genetic structure, phylogeographic patterns, diversity indices, and outlier detection. These codes were used to perform the analyses presented in the publication: https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.11407.



## SNP Discovery Workflow (GATK pipeline)

```mermaid
flowchart TD
    A["Use P. lilfordi reference genome (.fa)"] --> B["Create FASTA index and dictionary (Samtools)"]
    B --> C["Create BWA index of reference genome"]
    C --> D["Download FASTQs from Bassitta et al. 2021 via ENA Browser"]
    D --> E["Create sample sheet with sample IDs"]
    E --> F["Align FASTQ reads to reference genome (BWA)"]
    F --> G["Call variants with HaplotypeCaller (GATK)"]
    G --> H["Combine individual gVCFs (CombineGVCFs)"]
    H --> I["Filter variants (Get Supported Variants -D -N)"]
    I --> J["Generate VCF with GenotypeGVCFs (GATK)"]

## SNP Discovery Workflow (GATK pipeline) 
