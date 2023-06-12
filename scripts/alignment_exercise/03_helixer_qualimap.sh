#!/bin/bash 
#SBATCH --job-name=helixer_qualimap
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

##################################
# calculate stats on alignments
##################################
# this time we'll use qualimap



# input, output directories--------------------------------------------------------------

INDIR=../../results/alignment_exercise/alignments
OUTDIR=../../results/alignment_exercise/qualimap_reports/helixer
mkdir -p $OUTDIR

#get gff3 in gtf format
#module load singularity
#GFF=../../results/03_helixer/Arabidopsis_thaliana.gff3


#singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_convert_sp_gff2gtf.pl --gff $GFF -o ../../results/03_helixer/Arabidopsis_thaliana.gtf

# gtf annotation is required here
GTF=../../results/03_helixer/Arabidopsis_thaliana.gtf

# load software--------------------------------------------------------------------------
module load qualimap/2.2.1

qualimap \
    rnaseq \
    -bam $INDIR/SRR6852085.bam \
    -gtf $GTF \
    -outdir $OUTDIR/SRR6852085 \
    --java-mem-size=10G  

qualimap \
    rnaseq \
    -bam $INDIR/SRR6852086.bam \
    -gtf $GTF \
    -outdir $OUTDIR/SRR6852086 \
    --java-mem-size=10G  

##aggregate qualimap results
module load MultiQC/1.9

MultiQC_OUT=../../results/alignment_exercise/qualimap_reports/helixer/multiqc
multiqc -f -o $MultiQC_OUT $OUTDIR
