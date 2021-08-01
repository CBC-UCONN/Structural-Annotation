#!/bin/bash
#SBATCH --job-name=augutus
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

module load perl/5.32.1
export PERL5LIB=/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/3.01.03/PERL5.32/lib/perl5:/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/3.01.03/PERL5.32/lib/perl5/x86_64-linux-thread-multi:/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/3.01.03/PERL5.32/lib/perl5/x86_64-linux-thread-multi/auto:$PERL5LIB
export LD_LIBRARY_PATH=/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/3.01.03/PERL5.32/lib:$LD_LIBRARY_PATH

AUG_HOME=/labs/CBC/Tutorials/software/Augustus
export AUGUSTUS_CONFIG_PATH=${AUG_HOME}/config
export AUGUSTUS_BIN_PATH=${AUG_HOME}/bin
export PATH=${AUG_HOME}/bin:${AUG_HOME}/scripts:$PATH

MAKERDIR=/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/MAKER/3.01.03/bin
MAKERROUNDDIR=/labs/CBC/Tutorials/structural_annotation_for_assembled_genomes/13_marker/1_round/1_round_maker
MAKERROUND=first_iter
species=adb4.maker

if [[ -d $AUGUSTUS_CONFIG_PATH/species/$species ]]; then rm -r $AUGUSTUS_CONFIG_PATH/species/$species; fi


#take only the maker annotations
awk '{if ($2=="maker") print }' $MAKERROUNDDIR/$MAKERROUND.all.gff > maker_rnd1.gff


export AUGUSTUS_SCRIPTS=${AUG_HOME}/scripts

${AUGUSTUS_SCRIPTS}/gff2gbSmallDNA.pl maker_rnd1.gff ../1_round_maker/Athaliana_167_TAIR9.fa.masked 1000 $MAKERROUND.gb  # < elegans_genome_sm.fa 1000 $MAKERROUND.gb >

${AUGUSTUS_SCRIPTS}/randomSplit.pl $MAKERROUND.gb 100

${AUGUSTUS_SCRIPTS}/new_species.pl --species=$species

${AUGUSTUS_BIN_PATH}/etraining --species=$species --stopCodonExcludedFromCDS=true $MAKERROUND.gb.train

${AUGUSTUS_BIN_PATH}/augustus --species=$species $MAKERROUND.gb.test | tee firsttest.out


#evaluate the results 
#grep -A 22 Evaluation firsttest.out | awk {print $1} > eval_firsttest.out
## 


# optimize the model. this step could take a very long time
${AUGUSTUS_SCRIPTS}/optimize_augustus.pl --species=$species $MAKERROUND.gb.train --cpus=12 --kfold=12


#train again
${AUGUSTUS_BIN_PATH}/etraining --species=$species $MAKERROUND.gb.train
${AUGUSTUS_BIN_PATH}/augustus --species=$species $MAKERROUND.gb.test | tee optimizedtest.out

