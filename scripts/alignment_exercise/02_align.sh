#!/bin/bash
#SBATCH --job-name=hisat2_align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Align reads to genome
#################################################################
module load hisat2/2.2.1
module load samtools/1.12

INDIR=../../data/rnaseq/
OUTDIR=../../results/alignment_exercise/alignments
mkdir -p $OUTDIR

INDEX=../../results/alignment_exercise/hisat2_index/Athaliana


# run hisat2
hisat2 \
	-p 2 \
	-x $INDEX \
	-1 $INDIR/SRR6852085_1.fastq \
	-2 $INDIR/SRR6852085_2.fastq | \
samtools view -@ 1 -S -h -u - | \
samtools sort -@ 1 -T SRR6852085 - >$OUTDIR/SRR6852085.bam

# index bam files
samtools index $OUTDIR/SRR6852085.bam

# run hisat2 on second sample
hisat2 \
	-p 2 \
	-x $INDEX \
	-1 $INDIR/SRR6852086_1.fastq \
	-2 $INDIR/SRR6852086_2.fastq | \
samtools view -@ 1 -S -h -u - | \
samtools sort -@ 1 -T SRR6852086 - >$OUTDIR/SRR6852086.bam

# index bam files
samtools index $OUTDIR/SRR6852086.bam
