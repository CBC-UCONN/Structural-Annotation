#!/bin/bash
#SBATCH --job-name=hisat2_index
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

INDIR=../../results/repeatmasker/repeatmasker_out
GENOME=${INDIR}/GCF_000001735.4_TAIR10.1_genomic.fna.masked
INDEX=${INDIR}/Athaliana_masked

hisat2-build -p 8 ${GENOME} ${INDEX}

date


