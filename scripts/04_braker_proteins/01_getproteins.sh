#!/bin/bash
#SBATCH --job-name=get_proteins
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH --mem=10G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

OUTDIR=../../results/04_braker/proteins
    mkdir -p ${OUTDIR}
    cd ${OUTDIR}

# per the braker documentation, a protein database can be obtained as below
    # this code downloads protein sequences in fasta format for all genes in orthoDB for the Viridiplantae
    # then it concatenates them into a single fasta file
    # https://github.com/gatech-genemark/ProtHint#protein-database-preparation
wget --no-check-certificate https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz
tar -xvzf odb10_plants_fasta.tar.gz
cat plants/Rawdata/*fs >proteins.fa
