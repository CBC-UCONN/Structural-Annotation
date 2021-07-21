#!/bin/bash
#SBATCH --job-name=entap
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

##########################################
## EnTap				## 
##########################################
module load EnTAP/0.9.0-beta
module load diamond/0.9.36

EnTAP --runP \
	-i ../08_gFACs/mono_o/genes.fasta.faa \
	-d /isg/shared/databases/Diamond/RefSeq/complete.protein.faa.205.dmnd \
	-d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd \
	-d /isg/shared/databases/Diamond/ntnr/nr_protein.205.dmnd \
	--ontology 0 \
	--threads 8 

date


