#!/bin/bash
#SBATCH --job-name=04_repeatmasker
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH --mem=30G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

sed -i 's/Chr1.*/Chr1/g; s/Chr2.*/Chr2/g; s/Chr3.*/Chr3/g; s/Chr4.*4/Chr4/g; s/Chr5.*/Chr5/g; s/ChrM.*/ChrM/g; s/ChrC.*/ChrC/g;' Athaliana_167_TAIR9.fa

module load RepeatMasker/4.0.9-p2 
ln -s /labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/03_repeatmodeler/Athaliana_167_TAIR9.fa Athaliana_167_TAIR9.fa

#RepeatMasker -pa 8 -lib ../03_repeatmodeler/RM_150489.WedMar311224002021/consensi.fa -gff -a -noisy -xsmall Athaliana_167_TAIR9.fa
RepeatMasker -pa 8 -lib ../03_repeatmodeler/RM_*/consensi.fa -gff -a -noisy -xsmall Athaliana_167_TAIR9.fa
date
