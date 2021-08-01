#!/bin/bash
#SBATCH --job-name=maker_2
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

module load perl/5.32.1

script="/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/maker/bin/"

${script}/maker -base second_iter maker_opts.ctl maker_bopts.ctl maker_exe.ctl


