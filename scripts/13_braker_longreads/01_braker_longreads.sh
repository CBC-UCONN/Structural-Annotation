#!/bin/bash
#SBATCH --job-name=braker_longreads
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH --mem=20G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load BRAKER/2.1.6

# for testing new perl version:
module unload perl
module load perl/5.36.0

# set braker output directory and cd there
OUTDIR=../../results/braker/longreads
    mkdir -p ${OUTDIR}
cd ${OUTDIR}

# copy augustus config, set variable
export AUGUSTUS_CONFIG_PATH=$(pwd)/config
    cp -r /isg/shared/apps/augustus/3.6.0/config/ config

# create a temp directory for braker temp files
export TMPDIR=$(pwd)/braker_tmp 
    mkdir -p ${TMPDIR}

# run braker
BAM=../../alignments/mergednanopore.bam
GENOME=../../repeatmasker/repeatmasker_out/GCF_000001735.4_TAIR10.1_genomic.fna.masked

braker.pl --genome=${GENOME} \
        --bam ${BAM} \
        --softmasking 1 \
        --gff3 \
        --cores 16 \
        
date

