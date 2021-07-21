#!/bin/bash
#SBATCH --job-name=07_braker
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH --mem=20G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load python/3.6.3
module load biopython/1.70
module load N_GeneMark-ET/4.65_Perl5.32.1
module load samtools/1.10
module load bamtools/2.5.1
module load blast/2.10.0
module load genomethreader/1.7.1

BRAKER_HOME=/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/BRAKER/2.1.5
export PATH=${BRAKER_HOME}:${BRAKER_HOME}/scripts:${PATH}
AUG_HOME=/labs/CBC/Tutorials/software/Augustus
export AUGUSTUS_CONFIG_PATH=${AUG_HOME}/config
export AUGUSTUS_BIN_PATH=${AUG_HOME}/bin
export TMPDIR=$homedir/tmp
export BAMTOOLS_PATH=/isg/shared/apps/bamtools/2.5.1/bin
export BLAST_PATH=/isg/shared/apps/blast/ncbi-blast-2.10.0+/bin
export CDBTOOLS_PATH=/isg/shared/apps/qiime/1.9.1/cdbfasta/bin
export SAMTOOLS_PATH=/isg/shared/apps/samtools/1.9/bin
export DIAMOND_PATH=/isg/shared/apps/diamond/0.9.36/bin
#export PERL5LIB=/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/perl5/lib/site_perl/5.32.1/x86_64-linux-thread-multi:/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/perl5/lib/perl5:/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/BRAKER/2.1.5/scripts:$PERL5LIB


BAM=../06_short_read_align/finalbamfile.bam
GENOME=../04_repeatmasker/Athaliana_167_TAIR9.fa.masked

braker.pl --genome=${GENOME} \
        --bam ${BAM} \
        --softmasking 1 \
        --gff3 \
        --cores 16




date




