#!/bin/bash
#SBATCH --job-name=braker
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH --mem=20G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

# load modules, set environment variables to run braker
module load python/3.6.3
module load biopython/1.70
module load GeneMark-ET/4.68
module load samtools/1.10
module load bamtools/2.5.1
module load blast/2.10.0
module load genomethreader/1.7.1
module load ProtHint/2.6.0

# download braker from github
git clone https://github.com/Gaius-Augustus/BRAKER.git

BRAKER_HOME=BRAKER
export PATH=${BRAKER_HOME}:${BRAKER_HOME}/scripts:${PATH}
export AUGUSTUS_CONFIG_PATH=$(pwd)/config
    cp -r /isg/shared/apps/augustus/3.6.0_test/config/ config
export AUGUSTUS_BIN_PATH=/isg/shared/apps/augustus/3.6.0_test/bin/ 
export AUGUSTUS_SCRIPTS_PATH=/isg/shared/apps/augustus/3.6.0_test/scripts 
export TMPDIR=$(pwd)/braker_tmp 
    mkdir -p $TMPDIR
export BAMTOOLS_PATH=/isg/shared/apps/bamtools/2.5.1/bin
export BLAST_PATH=/isg/shared/apps/blast/ncbi-blast-2.10.0+/bin
export CDBTOOLS_PATH=/isg/shared/apps/qiime/1.9.1/cdbfasta/bin
export SAMTOOLS_PATH=/isg/shared/apps/samtools/1.9/bin
export DIAMOND_PATH=/isg/shared/apps/diamond/0.9.36/bin

# set variables for input/output files
GENOME=../04_repeatmasker/repeatmasker_out/GCF_000001735.4_TAIR10.1_genomic.fna.masked
PROTEINDB=proteins.fa

# run braker
braker.pl \
    --genome=${GENOME} \
    --prot_seq=${PROTEINDB} \
    --softmasking


