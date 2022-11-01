#!/bin/bash
#SBATCH --job-name=busco_filtered_split
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neranjan007@gmail.com
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              BUSCO                                   ##      
##########################################################

module load busco/5.0.0

AUG_HOME=/labs/CBC/Tutorials/software/Augustus
export AUGUSTUS_CONFIG_PATH=${AUG_HOME}/config

# output directories for busco results
OUTMONO=../../results/busco/braker_rna/monoexonics
    mkdir -p ${OUTMONO}

OUTMULTI=../../results/busco/braker_rna/multiexonics
    mkdir -p ${OUTMULTI}

# input directories for filtered protein sets
INMONO=../../results/gFACs/braker_rna/monoexonics
INMULTI=../../results/gFACs/braker_rna/multiexonics

# green plant busco database. already on the xanadu cluster. see documentation for how to obtain
BUSCODB=/isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10

# run busco
busco -i ${INMONO}/genes.fasta.faa \
        -o ${OUTMONO} \
        -c 8 \
        -l ${BUSCODB} \
        -m prot

busco -i ${INMULTI}/genes.fasta.faa \
        -o ${OUTMULTI} \
        -c 8 \
        -l ${BUSCODB} \
        -m prot


