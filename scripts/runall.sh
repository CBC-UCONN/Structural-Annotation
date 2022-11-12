#!/bin/bash

# get the data
cd 01_raw_data/
jid1=$( sbatch --parsable get_genome.sh )
jid2=$( sbatch --parsable sra_download.sh )

# repeat modeling/masking start when genome is downloaded
cd ../02_repeatmodeler
jid3=$( sbatch --parsable --dependency=afterok:$jid1 01_create_db.sh ) 
jid4=$( sbatch --parsable --dependency=afterok:$jid3 02_repeatmodeler.sh )

cd ../03_repeatmasker
jid5=$( sbatch --parsable --dependency=afterok:$jid4 01_repeatmasker.sh )
jid6=$( sbatch --parsable --dependency=afterok:$jid5 02_busco.sh )

# short read data processing, start when data are sequences are downloaded
cd ../04_fastq_qc
jid7=$( sbatch --parsable --dependency=afterok:$jid2 01_fastp.sh )

# short read alignment. start when repeatmasking and processing are done
cd ../05_short_read_align
jid8=$( sbatch --parsable --dependency=afterok:$jid7,$jid6 01_hisat2_index.sh )
jid9=$( sbatch --parsable --dependency=afterok:$jid8 02_align.sh )

# braker RNAseq, start when short reads are aligned
cd ../06_braker_RNAseq
jid10=$( sbatch --parsable --dependency=afterok:$jid9 01_braker_RNAseq )
jid11=$( sbatch --parsable --dependency=afterok:$jid10 02_busco.sh )

cd ../07_postprocess
jid12=$( sbatch --parsable --dependency=afterok:$jid10 01_gfacs_all.sh )
jid13=$( sbatch --parsable --dependency=afterok:$jid12 02_busco_all.sh )
jid14=$( sbatch --parsable --dependency=afterok:$jid12 03_EnTAP_all.sh )
jid15=$( sbatch --parsable --dependency=afterok:$jid12 04_EnTAP_monoexonics.sh )

# braker proteins, start when repeat masking is done
cd ../08_braker_proteins
jid16=$( sbatch --parsable --dependency=afterok:$jid6 01_getproteins.sh )
jid17=$( sbatch --parsable --dependency=afterok:$jid16 02_braker_proteins.sh )

cd ../09_postprocess
jid18=$( sbatch --parsable --dependency=afterok:$jid16 01_gfacs_all.sh )
jid19=$( sbatch --parsable --dependency=afterok:$jid18 02_busco_all.sh )
jid20=$( sbatch --parsable --dependency=afterok:$jid18 03_EnTAP_all.sh )
jid21=$( sbatch --parsable --dependency=afterok:$jid18 04_EnTAP_monoexonics.sh )

# TSEBRA, start when braker rnaseq and braker proteins are done
cd ../10_TSEBRA
jid22=$( sbatch --parsable --dependency=afterok:$jid10,$jid17 01_TSEBRA.sh )

cd ../11_postprocess
jid23=$( sbatch --parsable --dependency=afterok:$jid22 01_gfacs_all.sh )
jid24=$( sbatch --parsable --dependency=afterok:$jid23 02_busco_all.sh )
jid25=$( sbatch --parsable --dependency=afterok:$jid23 03_EnTAP_all.sh )
jid26=$( sbatch --parsable --dependency=afterok:$jid23 04_EnTAP_monoexonics.sh )

# braker longreads, start when data are downloaded and genome is masked
cd ../12_align_longreads
jid27=$( sbatch --parsable --dependency=afterok:$jid2,$jid6 01_minimap2.sh )

cd ../13_braker_longreads
jid28=$( sbatch --parsable --dependency=afterok:$jid27 01_braker_longreads.sh )
jid29=$( sbatch --parsable --dependency=afterok:$jid28 02_busco.sh )

cd ../14_postprocess
jid30=$( sbatch --parsable --dependency=afterok:$jid28 01_gfacs_all.sh )
jid31=$( sbatch --parsable --dependency=afterok:$jid30 02_busco_all.sh )
jid32=$( sbatch --parsable --dependency=afterok:$jid30 03_EnTAP_all.sh )
jid33=$( sbatch --parsable --dependency=afterok:$jid30 04_EnTAP_monoexonics.sh )

# gffcompare, start after all braker runs complete
cd ../15_gffcompare
jid34=$( sbatch --parsable --dependency=afterok:$jid10,$jid16,$jid22,$jid28 01_gffcompare.sh )

# orthofinder, start after all braker runs complete
cd ../16_orthofinder
jid35=$( sbatch --parsable --dependency=afterok:$jid10,$jid16,$jid22,$jid28 01_orthofinder.sh )
