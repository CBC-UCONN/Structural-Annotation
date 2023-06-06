#!/bin/bash
#SBATCH --job-name=EASEL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10G
#SBATCH --mail-user=

module load nextflow

SINGULARITY_TMPDIR=$PWD
export SINGULARITY_TMPDIR

nextflow run -hub gitlab PlantGenomicsLab/easel -profile singularity,xanadu -params-file params.yaml 


