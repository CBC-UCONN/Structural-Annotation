#!/bin/bash
#SBATCH --job-name=gffcmp
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load gffcompare/0.10.4

# output directory
OUTDIR=../../results/06_gffcompare
    mkdir -p ${OUTDIR}

# files
NCBIANNO=../../data/genome/GCF_000001735.4_TAIR10.1_genomic.gff
HELIXER=../../results/03_helixer/Arabidopsis_thaliana.gff3
BRAKERPROTEIN=../../results/04_braker/proteins/braker/augustus.hints.gtf
EASEL=../../results/05_EASEL/arabidopsis/final_predictions/arabidopsis_filtered.gtf

# run gffcompare
gffcompare -R -r <(awk '$7 !~ /?/' ${NCBIANNO}) -o ${OUTDIR}/helixer ${HELIXER}

gffcompare -R -r <(awk '$7 !~ /?/' ${NCBIANNO}) -o ${OUTDIR}/braker ${BRAKERPROTEIN}

gffcompare -R -r <(awk '$7 !~ /?/' ${NCBIANNO}) -o ${OUTDIR}/easel ${EASEL}

