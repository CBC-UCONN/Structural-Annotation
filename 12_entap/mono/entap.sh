#!/bin/bash
#SBATCH --job-name=entap
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neranjan007@gmail.com
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

##########################################
## EnTap                                ##
##########################################
module load EnTAP/0.9.0-beta
module load diamond/0.9.36

EnTAP --runP \
        -i ../../14_gfacs/2_round/mono_o/genes.fasta.faa \
        -d /isg/shared/databases/Diamond/RefSeq/complete.protein.faa.205.dmnd \
        -d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd \
        --ontology 0 \
        --threads 16

date

