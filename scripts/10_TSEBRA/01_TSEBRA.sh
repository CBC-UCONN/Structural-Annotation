#!/bin/bash
#SBATCH --job-name=TSEBRA
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=4G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load TSEBRA/1.0.3 

# output directory
OUTDIR=../../results/braker/tsebra
    mkdir -p ${OUTDIR}

# input files

# gene predictions
RNAGTF=../../results/braker/rnaseq/braker/augustus.hints.gtf
PROTGTF=../../results/braker/proteins/braker/augustus.hints.gtf

# hints files
RNAHINTS=../../results/braker/rnaseq/braker/hintsfile.gff
PROTHINTS=../../results/braker/proteins/braker/hintsfile.gff

# config
CFG=tsebra.cfg

# run TSEBRA
tsebra.py \
    -g ${RNAGTF},${PROTGTF} \
    -c ${CFG} \
    -e ${RNAHINTS},${PROTHINTS} \
    -o ${OUTDIR}/tsebra_combined.gtf

