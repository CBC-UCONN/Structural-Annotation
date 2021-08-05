#!/bin/bash
#SBATCH --job-name=gffcmp
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load gffcompare/0.10.4

gffcompare -R -r TAIR9_GFF3_genes.gff -o ref_braker braker_augustus.hints.gff3

gffcompare -R -r TAIR9_GFF3_genes.gff -o ref_marker marker_second_iter.all.gff

gffcompare -R -r TAIR9_GFF3_genes.gff -o ref_longRead longRead_braker.gtf


