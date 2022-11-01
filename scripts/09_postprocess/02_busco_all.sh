#!/bin/bash
#SBATCH --job-name=busco_all
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
date

##########################################################
##              BUSCO                                   ##      
##########################################################

module load busco/5.0.0

#export AUGUSTUS_CONFIG_PATH=$(pwd)/config
#    cp -r /isg/shared/apps/augustus/3.3.3/config/ config

BUSCODIR=../../results/busco/braker_proteins
mkdir -p ${BUSCODIR}
cd ${BUSCODIR}

OUTALL=all
OUTMONO=monoexonics
OUTMULTI=multiexonics

# green plant busco database. already on the xanadu cluster. see documentation for how to obtain
BUSCODB=/isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10

# predicted proteins
ALLPROT=../../../results/gFACs/braker_proteins/all/genes.fasta.faa
MONOPROT=../../../results/gFACs/braker_proteins/monoexonics/genes.fasta.faa
MULTIPROT=../../../results/gFACs/braker_proteins/multiexonics/genes.fasta.faa

# run busco
busco -i ${ALLPROT} \
        -o ${OUTALL} \
        -c 8 \
        -l ${BUSCODB} \
        -m prot \
        -f

busco -i ${MONOPROT} \
        -o ${OUTMONO} \
        -c 8 \
        -l ${BUSCODB} \
        -m prot \
        -f

busco -i ${MULTIPROT} \
        -o ${OUTMULTI} \
        -c 8 \
        -l ${BUSCODB} \
        -m prot \
        -f






