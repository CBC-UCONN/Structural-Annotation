#!/bin/bash
#SBATCH --job-name=getgenome
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=10G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

OUTDIR=../../data/genome
mkdir -p ${OUTDIR}

cd ${OUTDIR}

# Data: Arabidopsis thaliana TAIR10.1 assembly. 
    # NCBI Accession: GCF_000001735.4

# download genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz

# decompress
gunzip GCF_000001735.4_TAIR10.1_genomic.fna.gz

# strip off extra sequence name info after the space, e.g. 
    # this:         >NC_003070.9 Arabidopsis thaliana chromosome 1 sequence 
    # becomes this: >NC_003070.9
sed 's/ .*//' GCF_000001735.4_TAIR10.1_genomic.fna >tmp.fna
mv tmp.fna GCF_000001735.4_TAIR10.1_genomic.fna

# download the annotation (for comparison)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz

# decompress
gunzip GCF_000001735.4_TAIR10.1_genomic.gff.gz
