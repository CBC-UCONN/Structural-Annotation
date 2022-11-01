#!/bin/bash
#SBATCH --job-name=fastp_trimming
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Trimming/QC of reads using fastp
#################################################################
module load fastp/0.23.2

# set input/output directory variables
INDIR=../../data/rnaseq
OUTDIR=../../results/trimmed_fastqs
	mkdir -p ${OUTDIR}
REPORTDIR=../../results/trimmed_fastqs/fastp_reports
mkdir -p ${REPORTDIR}

# run fastp on accession SRR6852085
fastp \
	--in1 $INDIR/SRR6852085_1.fastq \
	--in2 $INDIR/SRR6852085_2.fastq \
	--out1 $OUTDIR/trim_SRR6852085_1.fastq \
	--out2 $OUTDIR/trim_SRR6852085_2.fastq \
	--json $REPORTDIR/SRR6852085_fastp.json \
	--html $REPORTDIR/SRR6852085_fastp.html

# run fastp on accession SRR6852086
fastp \
	--in1 $INDIR/SRR6852086_1.fastq \
	--in2 $INDIR/SRR6852086_2.fastq \
	--out1 $OUTDIR/trim_SRR6852086_1.fastq \
	--out2 $OUTDIR/trim_SRR6852086_2.fastq \
	--json $REPORTDIR/SRR6852086_fastp.json \
	--html $REPORTDIR/SRR6852086_fastp.html