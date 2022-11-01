!/bin/bash
#SBATCH --job-name=gfacs_general
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=40G
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


module load perl/5.32.1

genome="../04_repeatmasker/Athaliana_167_TAIR9.fa.masked"
alignment="../07_braker/braker/augustus.hints.gff3"
script="/labs/Wegrzyn/gFACs/gFACs.pl"

if [ ! -d general ]; then
        mkdir general
fi


perl "$script" \
        -f braker_2.1.2_gff3 \
        --statistics \
        --splice-table \
        --get-protein-fasta \
        --fasta "$genome" \
        -O general \
        "$alignment"


