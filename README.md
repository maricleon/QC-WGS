# QC-WGS
Pipeline for quality control of whole genome sequence data 

#Dependances: 

- minikraken2_v2_8GB_201904_UPDATE #database for kraken2
- NEBPE-PE.fa #adapter database if you are using NEBNext DNA prep kit for library preparation

#Note: please be aware that when utilizing the Illumina DNA preparation kit for library preparation, it is essential to ensure that the NEBPE-PE.fa file is appropriately adjusted to include the correct adapter sequences

#Programs:

- FastQC v0.11.9
- MultiQC v1.12
- TrimmomaticPE v0.39
- SPAdes v3.13.1
- Quast v5.2
- Bowtie2 v1.2.3
- Samtools 1.10
- Kraken2 v2.31 (minikraken2_v2 database)
- KronaTools 2.8.1


