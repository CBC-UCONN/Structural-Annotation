#!/bin/bash
#SBATCH --job-name=05_hisat2_index
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH --mem=20G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load hisat2/2.2.1

hisat2-build -p 8 ../04_repeatmasker/Athaliana_167_TAIR9.fa.masked Athaliana_masked


date


