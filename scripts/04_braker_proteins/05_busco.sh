#!/bin/bash
#SBATCH --job-name=busco_braker
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              BUSCO                                   ##      
##########################################################

module load busco/5.4.5

# green plant busco database. already on the xanadu cluster. see documentation for how to obtain
BUSCODB=/isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10

#Output directory
OUTDIR=../../results/04_braker/busco
mkdir -p ${OUTDIR}


# predicted proteins
PROTEINS=../../results/04_braker/agat/codingsequences/braker_proteins.fasta.faa


# run busco
busco -i ${PROTEINS} \
        -o ${OUTDIR} \
        -c 8 \
        -l ${BUSCODB} \
        -m prot \
        -f
