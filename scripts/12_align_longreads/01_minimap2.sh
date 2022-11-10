#!/bin/bash
#SBATCH --job-name=minimap2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=50G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

date
hostname

# load software
module load minimap2/2.24
module load samtools/1.16.1

# directories, files, etc
INDIR=../../data/nanopore
OUTDIR=../../results/alignments
mkdir -p ${OUTDIR}

PRESET="-ax splice -ub -k14"
GENOME=../../results/repeatmasker/repeatmasker_out/GCF_000001735.4_TAIR10.1_genomic.fna.masked

# run minimap2

# SRR10611193----------------------
SRR=SRR10611193
INFILE=${INDIR}/${SRR}.fastq
OUTFILE=${OUTDIR}/${SRR}.bam

minimap2 ${PRESET} ${GENOME} ${INFILE} -t 8 | \
    samtools view -@ 2 -S -h -u - | \
    samtools sort -@ 2 -T ${OUTDIR}/${SRR} - >${OUTFILE}

samtools index ${OUTFILE}

# SRR10611194----------------------
SRR=SRR10611194
INFILE=${INDIR}/${SRR}.fastq
OUTFILE=${OUTDIR}/${SRR}.bam

minimap2 ${PRESET} ${GENOME} ${INFILE} -t 8 | \
    samtools view -@ 2 -S -h -u - | \
    samtools sort -@ 2 -T ${OUTDIR}/${SRR} - >${OUTFILE}

samtools index ${OUTFILE}

# SRR10611195----------------------
SRR=SRR10611195
INFILE=${INDIR}/${SRR}.fastq
OUTFILE=${OUTDIR}/${SRR}.bam

minimap2 ${PRESET} ${GENOME} ${INFILE} -t 8 | \
    samtools view -@ 2 -S -h -u - | \
    samtools sort -@ 2 -T ${OUTDIR}/${SRR} - >${OUTFILE}

samtools index ${OUTFILE}


# merge bam files
samtools merge ${OUTDIR}/mergednanopore.bam ${OUTDIR}/SRR10611193.bam ${OUTDIR}/SRR10611194.bam ${OUTDIR}/SRR10611195.bam

# index bam file
samtools index ${OUTDIR}/mergednanopore.bam

# calculate stats
samtools stats ${OUTDIR}/mergednanopore.bam >${OUTDIR}/nanoporestats.txt