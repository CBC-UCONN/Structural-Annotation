#!/bin/bash
#SBATCH --job-name=entap_braker
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

##########################################
## EnTap                                ## 
##########################################
module load EnTAP/0.9.0-beta
module load diamond/0.9.36

OUTDIR=../../results/04_braker/EnTAP
    mkdir -p ${OUTDIR}
    
    cp entap_config.txt ${OUTDIR}
    cd ${OUTDIR}

PROTEINS=../agat/codingsequences/braker_proteins.fasta.faa

EnTAP --runP \
        -i ${PROTEINS} \
        -d /isg/shared/databases/Diamond/RefSeq/complete.protein.faa.216.dmnd \
        -d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd \
        --ontology 0 \
        --threads 16 

date
