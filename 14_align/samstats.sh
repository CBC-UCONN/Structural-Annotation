#!/bin/bash
#SBATCH --job-name=samstats
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mem=20G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

module load samtools/1.9

echo "SRR10611193"
samtools flagstat -@ 8 sorted_SRR10611193.bam
echo "SRR10611194"
samtools flagstat -@ 8 sorted_SRR10611194.bam
echo "SRR10611195"
samtools flagstat -@ 8 sorted_SRR10611195.bam

date


