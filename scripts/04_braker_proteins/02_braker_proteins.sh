#!/bin/bash
#SBATCH --job-name=braker_proteins
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH --mem=20G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load BRAKER/2.1.6
module load ProtHint/2.6.0

# for testing new perl version:
module unload perl
module load perl/5.36.0

# set braker output directory and cd there
OUTDIR=../../results/04_braker/proteins
    mkdir -p ${OUTDIR}
cd ${OUTDIR}

# copy augustus config, set variable
export AUGUSTUS_CONFIG_PATH=$(pwd)/config
    cp -r /isg/shared/apps/augustus/3.6.0/config/ config

# create a temp directory for braker temp files
export TMPDIR=$(pwd)/braker_tmp 
    mkdir -p ${TMPDIR}

# set variables for input/output files
GENOME="../../02_mask_repeats/repeatmasker/repeatmasker_out/GCF_000001735.4_TAIR10.1_genomic.fna.masked"
PROTEINDB=proteins.fa

# run braker
braker.pl \
    --genome=${GENOME} \
    --prot_seq=${PROTEINDB} \
    --softmasking \
    --cores 16 \
    --gff3
