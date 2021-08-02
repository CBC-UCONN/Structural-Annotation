#!/bin/bash
#SBATCH --job-name=gfacs_LR3
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=40G
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


module load perl/5.32.1

genome="../genome/Athaliana_167_TAIR9.fa.masked"
alignment="../02_braker/braker/augustus.hints.gff3"
script="/labs/Wegrzyn/gFACs/gFACs.pl"

if [ ! -d mono_o ]; then
        mkdir mono_o 
fi

if [ ! -d multi_o ]; then
        mkdir multi_o
fi

perl "$script" \
	-f braker_2.1.2_gff3 \
	--statistics \
        --statistics-at-every-step \
        --splice-table \
        --unique-genes-only \
        --rem-multiexonics \
        --rem-all-incompletes \
        --rem-genes-without-start-codon \
        --rem-genes-without-stop-codon \
	--min-CDS-size 300 \
        --get-protein-fasta \
        --fasta "$genome" \
        -O mono_o \
        "$alignment"

perl "$script" \
        -f braker_2.1.2_gff3 \
        --statistics \
        --statistics-at-every-step \
        --splice-table \
        --unique-genes-only \
        --rem-monoexonics \
        --min-exon-size 6 \
	--min-intron-size 9 \
	--min-CDS-size 300 \
	--get-protein-fasta \
        --fasta "$genome" \
        -O multi_o \
        "$alignment"


date


