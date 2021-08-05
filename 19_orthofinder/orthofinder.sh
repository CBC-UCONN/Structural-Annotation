#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH -c 16
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load OrthoFinder/2.5.1
module load muscle
module load DLCpar/1.0
module load FastME
module load diamond/0.9.25
module load mcl

orthofinder -f fasta -S diamond -t 16


