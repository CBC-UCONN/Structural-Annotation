#!/bin/bash
#SBATCH --job-name=gfacs_filtered_split
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=40G
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load gFACs/1.1.2

# output directories
OUTMONO=../../results/gFACs/braker_rna/monoexonics
mkdir -p $OUTMONO

OUTMULTI=../../results/gFACs/braker_rna/multiexonics
mkdir -p $OUTMULTI

# symlink genome
GENOME=../../results/repeatmasker/GCF_000001735.4_TAIR10.1_genomic.fna.masked
ln -s ${GENOME} genome.fa # symlink genome because gfacs does not support .fna suffix

# gff annotation
ANNOTATION=../../results/braker_rna/braker/augustus.hints.gff3

# calculate stats on monoexonics
gFACs.pl \
    -f braker_2.05_gff3 \
	--statistics \
	--statistics-at-every-step \
	--unique-genes-only \
	--rem-multiexonics \
	--rem-all-incompletes \
	--rem-genes-without-start-codon \
	--rem-genes-without-stop-codon \
	--get-protein-fasta \
	--fasta genome.fa \
	-O ${OUTMONO} \
	${ANNOTATION} 

# gFACs has a bug that causes it to occasionally print 2 fasta headers in a row for a sequence. 
# this line removes those extra headers
uniq $OUTMONO/genes.fasta.faa >genes.fasta.faa2 && mv genes.fasta.faa2 $OUTMONO/genes.fasta.faa

# calculate stats on multiexonics
gFACs.pl \
    -f braker_2.05_gff3 \
	--statistics \
	--statistics-at-every-step \
	--splice-table \
	--unique-genes-only \
	--rem-monoexonics \
	--rem-5prime-3prime-incompletes \
	--rem-genes-without-start-and-stop-codon \
	--get-protein-fasta \
	--min-exon-size 6 \
	--fasta genome.fa \
	-O ${OUTMULTI} \
	${ANNOTATION}

# gFACs has a bug that causes it to occasionally print 2 fasta headers in a row for a sequence. 
# this line removes those extra headers
uniq $OUTMULTI/genes.fasta.faa >genes.fasta.faa2 && mv genes.fasta.faa2 $OUTMULTI/general/genes.fasta.faa

date

