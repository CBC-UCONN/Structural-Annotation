#!/bin/bash
#SBATCH --job-name=busco
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              BUSCO                                   ##      
##########################################################

module load busco/5.4.5

# AUG_HOME=$HOME/Augustus
# export AUGUSTUS_CONFIG_PATH=${AUG_HOME}/config

# set output directory for busco, cd into parent directory
BUSCODIR=../../results/02_mask_repeats/busco/
mkdir -p ${BUSCODIR}
cd ${BUSCODIR}

OUTDIR=masked_genome

# input masked genome
MASKED_GENOME=../../02_mask_repeats/repeatmasker/repeatmasker_out/GCF_000001735.4_TAIR10.1_genomic.fna.masked

# green plant busco database, already downloaded on xanadu. see busco documentation for how to obtain
BUSCODB=/isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10

# run busco
busco -i ${MASKED_GENOME} \
        -o ${OUTDIR} \
        -c 8 \
        -f \
        -l ${BUSCODB} \
        -m genome


