#!/bin/bash

# get the data
cd 01_raw_data/
jid1=$( sbatch --parsable get_genome.sh )
jid2=$( sbatch --parsable sra_download.sh )

# repeat modeling/masking start when genome is downloaded
cd ../02_mask_repeats
jid3=$( sbatch --parsable --dependency=afterok:$jid1 01_create_db.sh ) 
jid4=$( sbatch --parsable --dependency=afterok:$jid3 02_repeatmodeler.sh )
jid5=$( sbatch --parsable --dependency=afterok:$jid4 03_repeatmasker.sh )
jid6=$( sbatch --parsable --dependency=afterok:$jid5 04_busco.sh )

#helixer
cd ../03_helixer
jid7=$( sbatch --parsable --dependency=afterok:$jid6 01_helixer.sh )
jid8=$( sbatch --parsable --dependency=afterok:$jid7 02_agat_stats.sh )
jid9=$( sbatch --parsable --dependency=afterok:$jid8 03_agat_codingsequences.sh )
jid10=$( sbatch --parsable --dependency=afterok:$jid9 04_busco.sh )
jid11=$( sbatch --parsable --dependency=afterok:$jid10 05_EnTAP.sh )

# braker
cd ../04_braker_proteins
jid12=$( sbatch --parsable --dependency=afterok:$jid6 01_getproteins.sh )
jid13=$( sbatch --parsable --dependency=afterok:$jid12 02_braker_proteins.sh )
jid14=$( sbatch --parsable --dependency=afterok:$jid13 03_agat_stats.sh )
jid15=$( sbatch --parsable --dependency=afterok:$jid14 04_agat_codingsequences.sh )
jid16=$( sbatch --parsable --dependency=afterok:$jid15 05_busco.sh )
jid17=$( sbatch --parsable --dependency=afterok:$jid16 06_EnTAP.sh )

# EASEL
cd ../05_EASEL
jid18=$( sbatch --parsable --dependency=afterok:$jid1,$jid2 easel.sh )

# compare
cd ../06_gffcompare
jid19=$( sbatch --parsable --dependency=afterok:$jid7,$jid13,$jid18 01_gffcompare.sh )
