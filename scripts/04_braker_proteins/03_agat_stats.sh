#!/bin/bash
#SBATCH --job-name=agat_stats_braker
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

OUTDIR="../../results/04_braker/agat/statistics"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sp_statistics.pl \
	--gff ../../proteins/braker/augustus.hints.gtf \
	-g ../../../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna \
	-o stats.txt


singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sq_stat_basic.pl \
	--gff ../../proteins/braker/augustus.hints.gtf \
	-g ../../../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna \
	-o basic_stats.txt

singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sp_functional_statistics.pl \
	--gff ../../proteins/braker/augustus.hints.gtf \
	-g ../../../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna \
	-o functional_statistics


