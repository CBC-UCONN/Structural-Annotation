#!/bin/bash
#SBATCH --job-name=repeatmodeler_model
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=50G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load RepeatModeler/2.0.4
module load ninja/0.95 

# go to repeat directory
REPDIR=../../results/02_mask_repeats

cd ${REPDIR}

# set repdb and run repeatmodeler
REPDB=athaliana_db

RepeatModeler -threads 30 -database ${REPDB} -LTRStruct 

date

