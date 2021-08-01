#!/bin/bash
#SBATCH --job-name=maker_1_analysis
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

#/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/maker/bin/maker -base first_iter maker_opts.ctl maker_bopts.ctl maker_exe.ctl


MAKERDIR="first_iter"
/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/maker/bin/maker -base ${MAKERDIR} -fix_nucleotides -dsindex

/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/maker/bin/gff3_merge -d ${MAKERDIR}.maker.output/${MAKERDIR}_master_datastore_index.log

/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/maker/bin/fasta_merge -d ${MAKERDIR}.maker.output/${MAKERDIR}_master_datastore_index.log

grep -c '>' ${MAKERDIR}.all.gff | awk {"print $1"} >  gene_annotations.txt

grep -c -P "\tgene\t" ${MAKERDIR}.all.gff | awk {"print $1"} > number_of_genes.txt

date

