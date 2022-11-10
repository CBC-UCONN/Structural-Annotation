#!/bin/bash
#SBATCH --job-name=maker_setup
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=4G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
