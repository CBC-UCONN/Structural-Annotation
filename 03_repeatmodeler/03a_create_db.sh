#!/bin/bash
#SBATCH --job-name=03_repeat_modeler_db
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH --mem=10G
#SBATCH -o log/%x_%j.out
#SBATCH -e log/%x_%j.err

hostname
date

module load RepeatModeler/2.01
ln -s /labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/03_repeatmodeler/Athaliana_167_TAIR9.fa Athaliana_167_TAIR9.fa

BuildDatabase -name "athaliana_db"  Athaliana_167_TAIR9.fa 
