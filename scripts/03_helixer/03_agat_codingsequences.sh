#!/bin/bash
#SBATCH --job-name=agat_codingsequences_helixer
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --mail-user=mia.nahom@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load singularity

OUTDIR="../../results/03_helixer/agat/codingsequences"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sp_keep_longest_isoform.pl \
	--gff ../../Arabidopsis_thaliana.gff3 \
	-o helixer_longest_isoform.gtf 


singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sp_extract_sequences.pl \
	-g helixer_longest_isoform.gtf \
	-f ../../../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna \
	-p \
	-o helixer_proteins.fasta.faa
