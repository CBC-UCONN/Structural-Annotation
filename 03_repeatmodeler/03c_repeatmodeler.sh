#!/bin/bash
#SBATCH --job-name=03b_repeatmodeler
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --partition=himem2
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH --mem=200G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load RepeatModeler/2.01
RepeatModeler -pa 30 -database athaliana_db -LTRStruct 

date

