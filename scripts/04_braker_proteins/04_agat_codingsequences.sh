#!/bin/bash
#SBATCH --job-name=agat_codingsequences_braker
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load singularity

OUTDIR="../../results/04_braker/agat/codingsequences"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sp_keep_longest_isoform.pl \
	--gff ../../proteins/braker/augustus.hints.gtf \
	-o braker_longest_isoform.gtf 


singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sp_extract_sequences.pl \
	-g braker_longest_isoform.gtf \
	-f ../../../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna \
	-p \
	-o braker_proteins.fasta.faa
