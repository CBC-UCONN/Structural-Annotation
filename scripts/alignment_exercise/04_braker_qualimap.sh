#!/bin/bash 
#SBATCH --job-name=braker_qualimap
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

##################################
# calculate stats on alignments
##################################
# this time we'll use qualimap

# load software--------------------------------------------------------------------------
module load qualimap/2.2.1


# input, output directories--------------------------------------------------------------

INDIR=../../results/alignment_exercise/alignments
OUTDIR=../../results/alignment_exercise/qualimap_reports/braker
mkdir -p $OUTDIR

# gtf annotation is required here
GFF=../../results/04_braker/proteins/braker/augustus.hints.gtf

qualimap \
    rnaseq \
    -bam $INDIR/SRR6852085.bam \
    -gtf $GFF \
    -outdir $OUTDIR/SRR6852085 \
    --java-mem-size=10G  

qualimap \
    rnaseq \
    -bam $INDIR/SRR6852086.bam \
    -gtf $GFF \
    -outdir $OUTDIR/SRR6852086 \
    --java-mem-size=10G  

##aggregate qualimap results
module load MultiQC/1.9

MultiQC_OUT=../../results/alignment_exercise/qualimap_reports/braker/multiqc
multiqc -f -o $MultiQC_OUT $OUTDIR
