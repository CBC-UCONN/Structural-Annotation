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

# output directory
OUTDIR=../../results/gffcompare
    mkdir -p ${OUTDIR}

# files
NCBIANNO=../../data/genome/GCF_000001735.4_TAIR10.1_genomic.gff
BRAKERRNA=../../results/braker/rnaseq/braker/augustus.hints.gff3
BRAKERPROTEIN=../../results/braker/proteins/braker/augustus.hints.gtf
BRAKERLONGREADS=../../results/braker/longreads/braker/augustus.hints.gff3

# run gffcompare
gffcompare -R -r <(awk '$7 !~ /?/' ${NCBIANNO}) -o ${OUTDIR}/rnaseq ${BRAKERRNA}

gffcompare -R -r <(awk '$7 !~ /?/' ${NCBIANNO}) -o ${OUTDIR}/proteins ${BRAKERPROTEIN}

gffcompare -R -r <(awk '$7 !~ /?/' ${NCBIANNO}) -o ${OUTDIR}/longreads ${BRAKERLONGREADS}


