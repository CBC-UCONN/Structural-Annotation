#!/bin/bash
#SBATCH --job-name=snap
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@gmail.com
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

module load perl/5.32.1

MAKERDIR=/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/maker/bin
MAKERROUNDDIR=/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/13_marker/1_round/1_round_maker

${MAKERDIR}/maker2zff ${MAKERROUNDDIR}/first_iter.all.gff

module unload perl/5.32.1

module load snap/2013-11-29
fathom -categorize 1000 genome.ann genome.dna

fathom -export 1000 -plus uni.ann uni.dna

#forge export.ann export.dna

hmm-assembler.pl first_iter . > first_iter.hmm

date


