#!/bin/bash
#SBATCH --job-name=align
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH --mem=50G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load hisat2/2.2.1
module load samtools/1.12

# set variables
INDIR=../../results/trimmed_fastqs
OUTDIR=../../results/alignments
	mkdir -p ${OUTDIR}

INDEX=../../results/repeatmasker/repeatmasker_out/Athaliana_masked

# run hisat2 for both sets of sequences
# set 1
SRR=SRR6852085
hisat2 \
	-p 8 \
	-x ${INDEX} \
	-1 ${INDIR}/trim_${SRR}_1.fastq \
	-2 ${INDIR}/trim_${SRR}_2.fastq | \
samtools view -@ 2 -S -h -u - | \
samtools sort -@ 2 -T ${OUTDIR}/${SRR} - >${OUTDIR}/${SRR}.bam

# index bam file
samtools index ${OUTDIR}/${SRR}.bam

# set 2
SRR=SRR6852086
hisat2 \
	-p 8 \
	-x ${INDEX} \
	-1 ${INDIR}/trim_${SRR}_1.fastq \
	-2 ${INDIR}/trim_${SRR}_2.fastq | \
samtools view -@ 2 -S -h -u - | \
samtools sort -@ 2 -T ${OUTDIR}/${SRR} - >${OUTDIR}/${SRR}.bam

# index bam file
samtools index ${OUTDIR}/${SRR}.bam

# merge bam files
samtools merge ${OUTDIR}/finalbamfile.bam ${OUTDIR}/SRR6852085.bam ${OUTDIR}/SRR6852086.bam

# index bam file
samtools index ${OUTDIR}/finalbamfile.bam

# calculate stats
samtools stats ${OUTDIR}/finalbamfile.bam >${OUTDIR}/stats.txt

date


