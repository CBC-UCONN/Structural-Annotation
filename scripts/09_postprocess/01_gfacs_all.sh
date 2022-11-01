#!/bin/bash
#SBATCH --job-name=gfacs_all
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=10G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load software
module load gFACs/1.1.2 

#output directory, variables
GFACS=../../results/gFACs
mkdir -p ${GFACS}

OUTALL=${GFACS}/braker_proteins/all
mkdir -p ${OUTALL}

OUTMONO=${GFACS}/braker_proteins/monoexonics
mkdir -p ${OUTMONO}

OUTMULTI=${GFACS}/braker_proteins/multiexonics
mkdir -p ${OUTMULTI}


# symlink genome
GENOME=../../results/repeatmasker/repeatmasker_out/GCF_000001735.4_TAIR10.1_genomic.fna.masked
ln -s ${GENOME} ${GFACS}/genome.fa # symlink genome because gfacs does not support .fna suffix

# gff annotation
#ANNOTATION=../../results/braker/proteins/braker/augustus.hints.gff3

# gtf annotation
ANNOTATION=../../results/braker/proteins/braker/augustus.hints.gtf

# run gfacs on all genes
gFACs.pl \
   -f braker_2.1.2_gtf \
   --statistics \
   --splice-table \
   --fasta ${GFACS}/genome.fa \
   --get-protein-fasta \
   -O ${OUTALL} \
   ${ANNOTATION}

# gFACs has a bug that causes it to occasionally print 2 fasta headers in a row for a sequence. 
# this line removes those extra headers
uniq ${OUTALL}/genes.fasta.faa >genes.fasta.faa2 && mv genes.fasta.faa2 ${OUTALL}/genes.fasta.faa


# calculate stats on monoexonics
gFACs.pl \
   -f braker_2.1.2_gtf \
	--statistics \
	--statistics-at-every-step \
	--unique-genes-only \
	--rem-multiexonics \
	--rem-all-incompletes \
	--rem-genes-without-start-codon \
	--rem-genes-without-stop-codon \
	--get-protein-fasta \
	--fasta ${GFACS}/genome.fa \
	-O ${OUTMONO} \
	${ANNOTATION} 

# gFACs has a bug that causes it to occasionally print 2 fasta headers in a row for a sequence. 
# this line removes those extra headers
uniq ${OUTMONO}/genes.fasta.faa >genes.fasta.faa2 && mv genes.fasta.faa2 ${OUTMONO}/genes.fasta.faa

# calculate stats on multiexonics
gFACs.pl \
   -f braker_2.1.2_gtf \
	--statistics \
	--statistics-at-every-step \
	--splice-table \
	--unique-genes-only \
	--rem-monoexonics \
	--rem-5prime-3prime-incompletes \
	--rem-genes-without-start-and-stop-codon \
	--get-protein-fasta \
	--min-exon-size 6 \
	--fasta ${GFACS}/genome.fa \
	-O ${OUTMULTI} \
	${ANNOTATION}

# gFACs has a bug that causes it to occasionally print 2 fasta headers in a row for a sequence. 
# this line removes those extra headers
uniq ${OUTMULTI}/genes.fasta.faa >genes.fasta.faa2 && mv genes.fasta.faa2 ${OUTMULTI}/genes.fasta.faa