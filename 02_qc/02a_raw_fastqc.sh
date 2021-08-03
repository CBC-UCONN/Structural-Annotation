#!/bin/bash
#SBATCH --job-name=raw_fastqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

##########################################################
## FASTQC Raw Reads                                     ##
##########################################################
module load fastqc/0.11.7
mkdir -p raw_fastqc_out

fastqc --threads 2 -o ./raw_fastqc_out ../01_raw_data/*.fastq 

