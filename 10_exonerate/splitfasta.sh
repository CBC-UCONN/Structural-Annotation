#!/bin/bash
#SBATCH --job-name=splitfasta
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=5G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load anaconda/4.4.0

echo `hostname`


protInPath=arabidopsis_files
protFileName=arabidopsis_filtered.pep

outdir="pieces"
mkdir -p pieces

script="splitfasta.py"

num_pieces=100

echo "Making pieces..."

python $script \
	--path $protInPath \
	--fasta $protFileName \
	--pathOut $outdir/ \
	--pieces $num_pieces

echo "Done!"

