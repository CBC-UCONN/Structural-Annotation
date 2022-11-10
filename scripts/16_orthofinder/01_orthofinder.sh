#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH -c 16
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load OrthoFinder/2.5.1
module load muscle
module load DLCpar/1.0
module load FastME
module load diamond/0.9.25
module load mcl

# create output directory 
OUTDIR=../../results/orthofinder
    mkdir -p ${OUTDIR}

# set up directory containing protein fasta, change suffix so orthofinder doesn't error out
FASTADIR=fasta
    mkdir -p ${OUTDIR}/${FASTADIR}

ln -s $(pwd)/../../results/braker/rnaseq/braker/augustus.hints.aa ${OUTDIR}/${FASTADIR}/augustus.hints.fa

# get arabidopsis lyrata proteins
wget -P ${OUTDIR}/${FASTADIR} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_protein.faa.gz
gunzip ${OUTDIR}/${FASTADIR}/GCF_000004255.2_v.1.0_protein.faa.gz

# run orthofinder
orthofinder -f ${OUTDIR}/${FASTADIR} -o ${OUTDIR}/braker_rnaseq -S diamond -t 16


