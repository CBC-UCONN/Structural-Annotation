#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Index the Genome
#################################################################

# load software
module load hisat2/2.2.0

# input/output directories
OUTDIR=../../results/alignment_exercise/hisat2_index
mkdir -p $OUTDIR

GENOME=../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna

hisat2-build -p 16 $GENOME $OUTDIR/Athaliana
