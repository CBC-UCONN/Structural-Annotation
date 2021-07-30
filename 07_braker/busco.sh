#!/bin/bash
#SBATCH --job-name=busco
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


busco -i ./braker/augustus.hints.aa \
        -o braker_aa \
        -c 8 \
        -l /isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10 -m prot


