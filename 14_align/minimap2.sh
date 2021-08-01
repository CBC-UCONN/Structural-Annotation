#!/bin/bash
#SBATCH --job-name=minimap2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=50G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


date
hostname

module load minimap2/2.17

preset="splice -uf -k14"
genome=../genome/Athaliana_167_TAIR9.fa.masked

minimap2 -ax $preset $genome ../data/SRR10611193_1.fasta -o SRR10611193.sam -t 8

minimap2 -ax $preset $genome ../data/SRR10611194_1.fasta -o SRR10611194.sam -t 8

minimap2 -ax $preset $genome ../data/SRR10611195_1.fasta -o SRR10611195.sam -t 8

