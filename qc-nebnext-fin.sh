#!/bin/bash

mkdir 1fastqc_pretrim
mkdir 2trimmomatic
mkdir 3fastqc_posttrim
mkdir 4spades_results
mkdir 5quast_results
mkdir 6kraken2
mkdir 7bowtie2
mkdir 8samtools

#dependances: minikraken2_v2_8GB_201904_UPDATE, NexteraPE-PE.fa or NEBPE-PE.fa, depends which library kit did we use


#first step is, to check the quality of raw reads before trimming
for file in ./*.fastq.gz; do
	fastqc $file ;
done

mv *.zip 1fastqc_pretrim;
mv *.html 1fastqc_pretrim;

#multiqc gatherd fastqc data into one report html or text files - your choice
~/.local/bin/multiqc ./1fastqc_pretrim/. -p #-p is for making a bunch of data in any case (if we analyse small amouint if data, it only generates html)
mv *.html 1fastqc_pretrim; 
mv multiqc_data 1fastqc_pretrim;
mv multiqc_plots 1fastqc_pretrim;

#now we trimm (we remove sequences shorter that 40 reads, adapters and bad quality sequences)
for file in $(ls ./*.fastq.gz | sed -r 's/_R[12]_001.fastq.gz//' | uniq); do
     TrimmomaticPE -threads 8 -phred33 "${file}_R1_001.fastq.gz" "${file}_R2_001.fastq.gz" "${file##*/}_R1_P1.trim.fastq.gz" "${file##*/}_R1_P2.trim.fastq.gz" "${file##*/}_R2_P1.trim.fastq.gz" "${file##*/}_R2_P2.trim.fastq.gz" ILLUMINACLIP:NEBPE-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:40;     
     echo ${file##*/}; 
done

rm *P2.trim.fastq.gz

#how about running fastqc now?
for file in ./*.trim.fastq.gz; do
	fastqc $file;
done

mv *.html 3fastqc_posttrim;
mv *.zip 3fastqc_posttrim;

#gathering fastqc files...
~/.local/bin/multiqc ./3fastqc_posttrim/. -p #-p da v vsakem primeru ven txt. ker če maš malo sekvenc jih ne da.
mv *.html 3fastqc_posttrim;
mv multiqc_data 3fastqc_posttrim;
mv multiqc_plots 3fastqc_posttrim;

#now, we embark on the assembly adventure!
for file in $(ls ./*.trim.fastq.gz | sed -r 's/_R[12]_P1.trim.fastq.gz//'| uniq);  do
	spades.py -1 "${file}_R1_P1.trim.fastq.gz" -2 "${file}_R2_P1.trim.fastq.gz" --careful --threads 8 --cov-cutoff auto -o "${file##/}_output"
	echo ${file##*/}
done

#this code enables us to rename an assembled FASTA file (fasta.contigs) according to the directory's name
for f in */contigs.fasta; do 
	cp $f ./${f%/contigs.fasta}.fasta;
	mv *.fasta 4spades_results;
done

#a little arrangement
mv *_output 4spades_results;
mv *.trim.fastq.gz 2trimmomatic;

#now, we assess the quality of the assembled sequences, including metrics such as the number of contigs, length, presence of Ns, GC content, N50, and more
~/quast-main/quast.py ./4spades_results/*.fasta
mv quast_results 5quast_results

eval "$(conda shell.bash hook)" #this is a trick to activate different conda environment inside the script
conda activate kraken #we use kraken for contamination check 

for file in $(ls ./*.fastq.gz | sed -r 's/_R[12]_001.fastq.gz//' | uniq); do
	kraken2 --threads 8 --db minikraken2_v2_8GB_201904_UPDATE --report "${file}_report" --gzip-compressed --paired "${file}_R1_001.fastq.gz" "${file}_R2_001.fastq.gz" > "${file}.kraken";
done

#to visualize kraken ouput (that comes as txt report) we can simply visualize it in more graphy way
#zazenem krona, aktiviram konda_env
#krona potebuje output krakena (ne report). z ukazom cut pa ohranimo samo 2 in 3 stolpec kar je edino potrebno za krono

eval "$(conda shell.bash hook)" #this is a trick to activate different conda environment inside the script 2
conda activate krona_env

for file1 in ./*.kraken; do
	cat $file1 | cut -f 2,3 > ${file1%.kraken}.krona; #with this command we use only 2nd and 3rd column from kraken report - its all we ever wanted
	ktImportTaxonomy ${file1%.kraken}.krona -o "${file1##*/}_html";
done

#a little arrangement
mv *.kraken 6kraken2;
mv *.krona 6kraken2;
mv *_html 6kraken2;
mv *_report 6kraken2;
rm -r *_html.files; #-r je brišeš mapo

#deativate all conda environments
conda deactivate;

#now we do bowtie - mapping reads onto assembled contigs

#first we index with bowtie2-build
for file in ./4spades_results/*_output.fasta; do
	bowtie2-build -f $file ${file%_output.fasta};
	mv ./4spades_results/*.bt2 7bowtie2;
	
done

#now, it's time for the real-deal action, not just the warm-up act
for file in $(ls ./2trimmomatic/*.trim.fastq.gz | sed -r 's/_R[12]_P1.trim.fastq.gz//' | uniq); do
	sample=`basename $file`
	bowtie2 --fr --end-to-end --threads 8 -x ./7bowtie2/${sample} 2>"${sample}_stats.txt" -sensitive -1 "${file}_R1_P1.trim.fastq.gz" -2 "${file}_R2_P1.trim.fastq.gz" -S "${sample}.sam";

done

mv *_stats.txt 7bowtie2
mv *.sam 7bowtie2


#samtools for estimating coverage
for file in ./7bowtie2/*.sam; do
	samtools view -S -b ${file} > "${file##*/}".bam;
	samtools sort "${file##*/}".bam -o "${file##*/}".sorted.bam;
	samtools depth "${file##*/}".sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${file##*/}.txt; 
done

#removing junk
rm ./7bowtie2/*.sam;
rm ./7bowtie2/*.bt2;
rm *.bam;
rm -r ./4spades_results/*_output;
mv *.txt 8samtools;
rm ./6kraken2/*krona;
rm ./6kraken2/*kraken


#this command edit each bowtie files in a way, that put in the first line, that is accualy name of a file. See where this is going
awk -i inplace -v ORS='\r\n' 'FNR==1{print FILENAME}1' ./7bowtie2/*_stats.txt;
#this command than gatherd all files in a final file, so the each two lines will contain name and value above
for FILE in ./7bowtie2/*.txt; do 
	sed -n '1p ; $p' $FILE > "${FILE}_temp";
done

cat ./7bowtie2/*_temp > ./7bowtie2/bowtie_fin.txt;
rm ./7bowtie2/*_temp;

#a simmilar for samtools. Naming first line after files name and than gathers all data in final file
awk -i inplace -v ORS='\r\n' 'FNR==1{print FILENAME}1' ./8samtools/*.sam.txt
cat ./8samtools/*.txt > samtools_fin
mv samtools_fin 8samtools
