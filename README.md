# QC-WGS
Pipeline for quality control of whole genome sequence data 

#Dependances: 

- minikraken2_v2_8GB_201904_UPDATE #database for kraken2
- NexteraPE-PE.fa #adapter database if you are using Illumina DNA prep kit for library preparation 
- NEBPE-PE.fa #adapter database if you are using NEBNext DNA prep kit for library preparation

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

#after analysis, from platform directory run:

sed -e n\;d <./1fastqc_pretrim/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt > stevilo_sekvenc_pretrim
mv stevilo_sekvenc_pretrim ./1fastqc_pretrim/multiqc_data/

sed -e n\;d <./3fastqc_posttrim/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt > stevilo_sekvenc_posttrim
mv stevilo_sekvenc_posttrim ./3fastqc_posttrim/multiqc_data/

