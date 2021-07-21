#!/bin/bash
#SBATCH --job-name=exonerate
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=1-100%20
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o ./logfiles/%A_%a_%j.out
#SBATCH -e ./logfiles/%A_%a_%j.err

hostname
date

##########################################
## 				## 
##########################################
module load exonerate/2.4.0

genome=../04_repeatmasker/Athaliana_167_TAIR9.fa.masked
species="arabidopsis"

mkdir -p gffpieces
out=gffpieces/${species}_${SLURM_ARRAY_TASK_ID}.exon.gff
pep=pieces/arabidopsis_filtered.pep${SLURM_ARRAY_TASK_ID}.fa


exonerate --model protein2genome \
	--query $pep \
	--target $genome \
	-n 1 --percent 95 --score 500 --minintron 50 --showalignment no \
	--showtargetgff yes --geneseed 250 --forcegtag \
	--hspfilter 100 --showvulgar yes --maxintron 200000 > $out

