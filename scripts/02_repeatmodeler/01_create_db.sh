#!/bin/bash
#SBATCH --job-name=repeatmodeler_db
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=10G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load RepeatModeler/2.01

# set variables
REPDIR=../../results/repeatmodeler
    mkdir -p ${REPDIR}

# genome
GENOME=../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna 

BuildDatabase -name "${REPDIR}/athaliana_db" ${GENOME}
