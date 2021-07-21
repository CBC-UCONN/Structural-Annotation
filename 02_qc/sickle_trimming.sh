#!/bin/bash
#SBATCH --job-name=sickle_trimming
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


module load sickle/1.33

sickle pe \
        -t sanger \
        -f ../01_raw_data/SRR6852085_1.fastq -r ../01_raw_data/SRR6852085_2.fastq \
        -o trimmed_SRR6852085_1.fastq -p trimmed_SRR6852085_2.fastq -s trimmed_singles_6852085.fastq \
        -q 30 -l 50


sickle pe \
        -t sanger \
        -f ../01_raw_data/SRR6852086_1.fastq -r ../01_raw_data/SRR6852086_2.fastq \
        -o trimmed_SRR6852086_1.fastq -p trimmed_SRR6852086_2.fastq -s trimmed_singles_6852086.fastq \
        -q 30 -l 50 


date
