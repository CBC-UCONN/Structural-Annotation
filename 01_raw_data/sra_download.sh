#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=20G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load sratoolkit

fastq-dump --split-files SRR6852085
fastq-dump --split-files SRR6852086

date
