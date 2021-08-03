#!/bin/bash
#SBATCH --job-name=busco
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
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

AUG_HOME=/labs/CBC/Tutorials/software/Augustus
export AUGUSTUS_CONFIG_PATH=${AUG_HOME}/config

busco -i general/genes.fasta.faa \
        -o 2_round_maker_general \
        -c 8 \
        -l /isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10 -m prot

busco -i mono_o/genes.fasta.faa \
        -o 2_round_maker_mono \
        -c 8 \
        -l /isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10 -m prot

busco -i multi_o/genes.fasta.faa \
        -o 2_round_maker_multi \
        -c 8 \
        -l /isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10 -m prot


