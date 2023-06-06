#!/bin/bash
#SBATCH --job-name=repeatmasker
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

# load software
module load RepeatMasker/4.1.2

# set variables
REPLIB=../../results/02_mask_repeats/athaliana_db-families.fa
OUTDIR=../../results/02_mask_repeats/repeatmasker
    mkdir -p ${OUTDIR}

GENOME=../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna

RepeatMasker -dir repeatmasker_out -pa 8 -lib ${REPLIB} -gff -a -noisy -xsmall ${GENOME}

date
