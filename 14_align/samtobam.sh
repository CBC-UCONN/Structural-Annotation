#!/bin/bash
#SBATCH --job-name=samtobam
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=20G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

module load samtools/1.9

prefix="arabidopsis"

for file in *.sam
do
        fname=${file#../sam/}
        echo $fname
        samtools view -b -@ 16 $file | samtools sort -o sorted_${fname%.sam}.bam -@ 16

done


samtools merge ${prefix}.merged.bam sorted_*.bam


