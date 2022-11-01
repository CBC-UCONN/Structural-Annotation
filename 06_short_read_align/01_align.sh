#!/bin/bash
#SBATCH --job-name=06_align
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH --mem=50G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load hisat2/2.2.1

hisat2 -x ../05_hisat2_index/Athaliana_masked \
	-1 ../02_qc/trimmed_SRR6852085_1.fastq -2 ../02_qc/trimmed_SRR6852085_2.fastq \
	-p 8 \
	-S SRR6852085.sam 


hisat2 -x ../05_hisat2_index/Athaliana_masked \
	-1 ../02_qc/trimmed_SRR6852086_1.fastq -2 ../02_qc/trimmed_SRR6852086_2.fastq \
	-p 8 \
	-S SRR6852086.sam 


module unload hisat2/2.2.1

module load samtools/1.10 

samtools view -@ 8 -uhS SRR6852085.sam | samtools sort -@ 8 -T SRR6852085 -o sorted_SRR6852085.bam
samtools view -@ 8 -uhS SRR6852086.sam | samtools sort -@ 8 -T SRR6852086 -o sorted_SRR6852086.bam  


samtools merge finalbamfile.bam sorted_SRR6852085.bam sorted_SRR6852086.bam

echo "sorted_SRR6852085.bam"
samtools flagstat -@ 8 sorted_SRR6852085.bam 

echo "sorted_SRR6852086.bam"
samtools flagstat -@ 8 sorted_SRR6852086.bam

date


