#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=8G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# run this script to download the data from the SRA
    # if you're on xanadu and want to save time, run the other
    # script in this directory to symlink the data instead

# data
# Illumina paired-end RNA-seq from Arabidopsis leaf tissue
    # bioproject: PRJNA438701
        # biosample: SAMN08724106
        # SRA runs:
            # SRR6852085
            # SRR6852086

# Oxford nanopore long read cDNA
    # bioproject: PRJNA594286 
        # biosamples: SAMN13510394,SAMN13510392,SAMN13510391
        # SRA runs:
            # SRR10611193
            # SRR10611194
            # SRR10611195

# load software
module load sratoolkit/2.11.3

# output directory, create if it doesn't exist
OUTDIR=../../data/rnaseq
mkdir -p ${OUTDIR}

cd ${OUTDIR}

fasterq-dump SRR6852085
fasterq-dump SRR6852086

# nanopore outdir
NANODIR=../nanopore
mkdir -p ${NANODIR}

cd ${NANODIR}

fasterq-dump SRR10611193
fasterq-dump SRR10611194
fasterq-dump SRR10611195

date
