#!/bin/bash
#SBATCH --job-name=symlink_fastqs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=8G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# run this script to symlink the data instead of downloading from the SRA
    # if you want to download the data instead, run the other script

# data
# RNA-seq from Arabidopsis leaf tissue
    # bioproject: PRJNA438701
        # biosample: SAMN08724106
        # SRA runs:
            # SRR6852085
            # SRR6852086

# output directory, create if it doesn't exist
OUTDIR=../../data/rnaseq
mkdir -p $OUTDIR

cd ${OUTDIR}

# symlink data from fpath
fpath="/core/cbc/tutorials/rawdata/Structural_Annotation/v1/arabidopsis_leaf/"

for f in ${fpath}*; do
        echo $f
        echo `basename ${f}`
        ln -s ${f} `basename ${f}`
done



# Oxford nanopore long read cDNA
    # bioproject: PRJNA594286 
        # biosamples: SAMN13510394,SAMN13510392,SAMN13510391
        # SRA runs:
            # SRR10611193
            # SRR10611194
            # SRR10611195

# output directory, create if it doesn't exist
OUTDIR=../../data/nanopore
mkdir -p $OUTDIR

cd ${OUTDIR}

# symlink data from fpath
fpath="/core/cbc/tutorials/rawdata/Structural_Annotation/v1/nanopore/"

for f in ${fpath}*; do
        echo $f
        echo `basename ${f}`
        ln -s ${f} `basename ${f}`
done
