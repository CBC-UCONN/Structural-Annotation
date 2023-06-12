This repository is a tutorial for genome annotation. The scripts are set up to run on UConn's Xanadu cluster, including Xanadu-specific SLURM headers and software modules. To run it on Xanadu, simply clone this repository and start submitting the scripts by following along with this readme. If you are interested in running it elsewhere, you'll need to install the relevant software and alter or remove the SLURM headers, but otherwise, the tutorial is self-contained and pulls the necessary data from public databases.

Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use `sbatch scriptname` after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as `nano`, `vim`, or `emacs`.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  The tutorial assumes basic familiarity with Linux. In this guide, you will be working with common bioinformatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.

Contents
1.   [Overview](#1-overview)
2.   [Downloading Data](#2-downloading-the-data)
3.   [Identifying and Masking Repetitive Elements](#3-identifying-and-masking-repetitive-elements)
4.   [Helixer](#4-helixer-annotation) 
5.   [Helixer Evaluation](#5-helixer-evaluation)  
6.   [Braker Annotation with Protein Evidence](#6-braker-annotation-with-protein-evidence)  
7.   [Braker Evaluation](#7-braker-evaluation)   
8.   [EASEL Pipeline](#8-easel-pipeline)  
9.   [Compare Evaluations with gffcompare](#9-compare-evaluations)   

## 1.  Overview
In this tutorial, we'll annotate an *Arabidopsis thaliana* genome using a variety of methods. We'll use both RNA-seq data from leaf tissue and a database of green plant proteins as evidence in the annotation process. To get started, you can clone the tutorial repository to obtain all the scripts we're going to run. If you're working on UConn's Xanadu cluster, you can submit each script without modification. If you're working elsewhere, you will have to make sure the software is installed and modify or remove the SLURM headers, as appropriate to your system. To clone the repository:

```
git clone git@github.com:CBC-UCONN/Structural-Annotation.git
```

You'll notice the structure of the directory is initially pretty spare. Each script will create, as necessary, directories to hold data and results as you move along. 

### `SLURM`
UConn's Xanadu cluster uses `SLURM` to manage resources. Scripts to accomplish some task or requests for resources for interactive analysis, are routed through it. If you're working through the tutorial on Xanadu, you can submit each script using the command `sbatch` (e.g. `sbatch get_genome.sh`). `SLURM` will then queue your job, and when resources become available, run your code. Any output from your code that is written to the *standard output* or *standard error* channels is captured for later review. When submitting scripts via `sbatch`, we usually specify the resource request in a *header*. We'll review an example header here so we don't have to for every script going forward. 

```bash
#!/bin/bash
#SBATCH --job-name=getgenome
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@uconn.edu
```

The first line is the *shebang*, indicating the script is `bash` code. The following lines, beginning with `#SBATCH`, pertain to `SLURM`. We specify a job name with `--job-name`, and the format of the names of the files that capture the standard output (`-o`) and the standard error (`-e`). In this case for `-o`, the format will be `JOBNAME_JOBID.out`. You should always check the `.out` and `.err` files after a job completes. These will usually contain any error messages produced by your code. 

The next six lines specify the resources requested. `-N 1` asks that the job be run on 1 node. `-n 1` tells SLURM we're only going to launch 1 task. `-c 1` requests 1 CPU for this task. `--mem=10G` requests 10 gigabytes of memory for this task. `--partition` asks that this job be run on the general partition (we also have high memory and GPU partitions) and `--qos=general` specifies the "general" quality of service. For most applications, we only vary `-c`, `--mem`, `--partition`, and `--qos=general`. 

It's important that you request adequate, but not excessive resources for each task. Also, generally you must tell the software you are running how many CPUs it can use (and sometimes how much memory) or it may not use all that you have requested. You can read more about resource allocations on Xanadu [here](https://github.com/CBC-UCONN/CBC_Docs/wiki/Requesting-resource-allocations-in-SLURM) and at the link above. 

After the SLURM header, we have the code to be run, as you will see in the sections below. 

### Software modules

You will note that in most of these scripts right after the SLURM header there is a `module load` command (e.g. `module load sratoolkit/2.11.3`). On Xanadu, we use a module system that allows us to have lots of software installed, often with conflicting dependencies, and to maintain functional older versions of software. In order to use one of these pieces of software, however, it must first be loaded. If you are working through this code on another system, you may need to install the software yourself, or modify the way it is made available (e.g. the version numbers). 

##2 Downloading the Data

The first step in the tutorial is to obtain the data we need. There are three main pieces of data: a reference genome, RNA-seq data, and protein data from green plants. We'll get the genome and RNA-seq data now, and download the protein data when we use it later. 

Scripts to obtain the data can be found in the directory [`scripts/01_raw_data`](). There are three. One downloads our reference genome, the other two are approaches for getting the RNA-seq data. 

### Downloading the genome

To get the genome, enter that directory and run the script `get_genome.sh`. By typing `sbatch get_genome.sh` on the command line. The body of the script (after the SLURM header) looks like this: 

```bash
OUTDIR=../../data/genome
mkdir -p ${OUTDIR}

cd ${OUTDIR}

# Data: Arabidopsis thaliana TAIR10.1 assembly. 
    # NCBI Accession: GCF_000001735.4

# download genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz

# decompress
gunzip GCF_000001735.4_TAIR10.1_genomic.fna.gz

# strip off extra sequence name info after the space, e.g. 
    # this:         >NC_003070.9 Arabidopsis thaliana chromosome 1 sequence 
    # becomes this: >NC_003070.9
sed 's/ .*//' GCF_000001735.4_TAIR10.1_genomic.fna >tmp.fna
mv tmp.fna GCF_000001735.4_TAIR10.1_genomic.fna
```

We first create a variable and assign it the name of our desired output directory. We then create that directory and `cd` into it. Most scripts in this tutorial will begin by specifying variables for input and output directories (and sometimes other relevant files and directories). Not all will move into the output directory for the analysis, but direct output there instead. 

The next steps are to download the genome from NCBI, decompress it and edit the sequence names, which have extraneous information (everything after the first space) that can trip up some software. You will note after this script completes, that there is now a new directory `data/genome` in the root of the tutorial repository containing our genome file (in this case with the suffix `.fna`). 

The reference genome is in *FASTA* format. This format is plain text. Each sequence in the file begins with a header line, followed by the sequence:

```
>NC_003070.9
ccctaaaccctaaaccctaaaccctaaacctctGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCATG
AATCCCTAAATACCTAAttccctaaacccgaaaccggTTTCTCTGGTTGAAAATCATTGTGtatataatgataattttat
CGTTTTTATGTAATTGCTTATTGTTGTGtgtagattttttaaaaatatcatttgagGTCAATACAAATCCTATTTCTTGT
GGTTTTCTTTCCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATT
...
```

### The RNA-seq data

There are two scripts that can be used to obtain the RNA-seq data. If you are working on Xanadu, you can use `symlink_data.sh`. That script will create a new directory, `data/rnaseq` and create *symlinks*, or pointers that point to copies of the data already residing on the cluster. If you would rather download the data from NCBI, you can run `sra_download.sh`. The body of the script looks like this:

```bash
# RNA-seq from Arabidopsis leaf tissue
    # bioproject: PRJNA438701
        # biosample: SAMN08724106
        # SRA runs:
            # SRR6852085
            # SRR6852086

# load software
module load sratoolkit/2.11.3

# output directory, create if it doesn't exist
OUTDIR=../../data/rnaseq
mkdir -p $OUTDIR

cd ${OUTDIR}

fasterq-dump SRR6852085
fasterq-dump SRR6852086
```

Here we are using `sratoolkit` and the command `fasterq-dump` to download two different sets of RNA-seq data derived from Arabidopsis leaf tissue. 

The resulting data will be found in the newly created `data/rnaseq` directory. There will be four files of sequence data, each in *fastq* format:

```
SRR6852085_1.fastq
SRR6852085_2.fastq
SRR6852086_1.fastq
SRR6852086_2.fastq
```

Our data are *paired-end*, which means that each molecule of cDNA in the sequence library is sequenced twice, once from each end. So our sequence data comes in pairs of files (`_1.fastq` and `_2.fastq`). The mate-pairs are kept in order, so if you do anything to disrupt the sequence ordering of these files that information will be scrambled. 

*fastq* format is a bit like *fasta*, except that it contains more information. Each sequence record has four lines: a header line giving the sequence name (beginning with `@`), a nucleotide sequence, a comment line beginning with `+`, and a line containing ASCII encoded, phred-scaled base qualities. For more on base qualities see [here](https://www.drive5.com/usearch/manual/quality_score.html). 

```
@SRR6852085.1 1 length=150
NGGCATGCAGACTTGTGAGGGGACGGAGACAATTCCGCTCANNGCNAGGTCACACACGTGTCTATTGTNAGGTGTGTACATNGGCAACGTGAAAGCGATAGTGAGGGNACAGTTTGGAATGTACAGCTGAACGGACATCACACGAAACCT
+SRR6852085.1 1 length=150
#<AFFFJJJJJJJF-FF7FFAJAJJJJJJJJFFJJJJFFJ7##FJ#<JJFJJJJJJJA-AAJFAAAAA#AJF7F-FAJJJJ#JJ<AFJJA-<F<----7----7-J-#-777-7AAF-AAF-FAJ--7--7--7--7-7--7--------
```

You can run this script by entering the directory `scripts/01_raw_data/` and typing `sbatch sra_download.sh` on the command line. 

Once the data is downloaded the folder will have the following files:
```
01_raw_data/
├── SRR6852085_1.fastq
├── SRR6852085_2.fastq
├── SRR6852086_1.fastq
└── SRR6852086_2.fastq
```

You can check the number of reads in a fastq file using the following awk command:
```
awk '{s++}END{print s/4}' SRR6852085_1.fastq
```
As each fastq file contains 4 lines per each read, you will need to divide the total number of counts you get, by that number, which will give you the number of reads in the sample.   

Additionally, there is a script labeled `symlink_data.sh`. If you are working on UConn's Xanadu cluster you can run this script to access the data, rather than downloading a copy from NCBI.

## 3. Indentifying and Masking Repetitive Elements

The first step in annotating a genome is to identify repetitive elements. While this is a worthy goal in its own right, repetitive elements can also interfere with gene prediction, so by finding them first, we can mask them out and do a much better job finding genes. The process we are going to go through here can be thought of as a first pass at finding repetitive elements. The repeat models we generate will be ok, but highly redundant, and many will be fragmented. Getting a good candidate set of repetitive elements and understanding the distribution of each type throughout the genome requires more manual curation than we will undertake here. 

There are two steps in this phase, first we're going to discover repetitive elements in our genome using `RepeatModeler`, and then we will use the resulting library of repeats to mask the genome. We'll use the masked genome for gene prediction in the next steps. 

### Discovering repetitive elements with `RepeatModeler`

`RepeatModeler` is a pipeline that ties together a few different pieces of software, but because repeats aren't a major focus of this analysis, we won't go into detail here. The main thrust is that we'll be identifying repetitive elements *de novo*, i.e. without using a library of elements from a related species. 

A significant proportion of any eukaryotic genome is low complexity regions, often consisting of [repetitive elements](https://en.wikipedia.org/wiki/Repeated_sequence_(DNA)). Identifying these repeats is an important part of structural annotation of genomic sequence. A good fraction of repeats are associated with Transposable elements (TE) also known as mobile elements. TE are biologically important as they are proposed to have role in [genome evolution](https://pubmed.ncbi.nlm.nih.gov/15016989/) , [genomic rearrangement](https://pubmed.ncbi.nlm.nih.gov/15020798/) and modulation of gene expression. Repeats can negatively impact gene prediction, so we need to identify, annotate, and mask them.  There are two ways to mask sequence:


(1) **Hard-Masking** where the the sequence is replaced by Ns.  Example the repeat sequence TGCAAATCGCA (terminal inverted repeat sequence  of Class 2 TE's) is hard masked in the sequence below

`CTGTGCAAATCGCAGTTA -> CTGNNNNNNNNNNNGTTA `

(2) **Soft-masking** , this involves converting the sequence from Uppercase to lowercase as an examle the same repeat sequence is soft masked below

`CTGTGCAAATCGCAGTTA -> CTGtgcaaatcgcagTTA `

For our downstream software, softmasking is preferable. It will let annotators know that a given sequence is repetitive, so that it will be ignored in initial gene-finding, but retaining the sequence allows nearby genic sequence to be extended into it if necessary.

Getting good annotations of repetitive elements requires a fair bit of manual curation of repeat element models. 

Software used in identification of repeats can be categorised as extrensic and intrinsic tools.

**Extrinsic tools**, e.g. [RepeatMasker](https://www.repeatmasker.org/), uses the repeat sequence (from a closely related species) listed in Repbase (a repeat database) and annotate there presence in our assembled genome.

**Intrinsic tool**, e.g. [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/), perform a de novo identification and modelling of TE families. 

### RepeatModeler
It relies on three de-novo repeat finding programs ( RECON, RepeatScout and LTRHarvest/LTR_retriever ) which employ complementary computational methods for identifying repeat element boundaries and family relationships from sequence data. A brief description each of the programme is given below

**RECON** can carry out denovo identification and classification of repeat sequence families. In order to achieve that it first perform pairwise alignment between the genomic sequences and then it cluster the sequences using single linkage cluster (agglomerative clustering) approach. It then uses the multiple alignment information to define boundaries of repeats and also to distinguish homologous but distinct repeat element families.

**RepeatScout** first create a frequency table of l-mers (l=ceil(log_4(L)+1), where L length of input sequence) then it creates a fasta file of all the repetitive elements. It then run 2 rounds of filtering on the repeat elements ("filter-stage-1.prl" and "filter-stage-2.prl"), in the first round it removes low complexity tandem repeats from the fasta file following that the frequency of repeats in fasta fasta file is estimated in the genome using RepeatMasker. In second round of filtering, repeats appearing less than certain number of times (default 10) are removed.

**LTRHarvest/LTR_retriever** perform de novo detection of full length LTR retrotransposons in large sequence sets. LTRharvest efficiently delivers high quality annotations based on known LTR transposon features like length, distance, and sequence motifs. LTR_retriever carry out accurate identification of LTR retrotransposons (LTR-RTs) from output of LTRharvest and generates non-redundant LTR-RT library for genome annotations.

Typically both Intrinsic and Extrinsic tools are used to annotate and mask the repaets of a genome and thats what we will be doing here. We will perform Repeatmodeller first to identify novel motifs and then using RepeatMasker we will softmask repeats in the genome using sequences of novel repeats and reference repeats (from Repbase).  Before we identify our repeat regions, we must first compile our database using the "BuildDatabase" command of RepeatModeler. This will format the FASTA files for use with RepeatModeler.
```
module load RepeatModeler/2.0.4
BuildDatabase -name "athaliana_db"  Athaliana_167_TAIR9.fa
```
Command options:
```
BuildDatabase [-options] -name "mydb" <seqfile(s) in fasta format>
-name <database name>  The name of the database to create.
```
The complete slurm script called 01_create_db.sh can be found in `02_mask_repeats/` directory.
This will create the following database files:
```
├── athaliana_db.nhr
├── athaliana_db.nin
├── athaliana_db.nnd
├── athaliana_db.nni
├── athaliana_db.nog
├── athaliana_db.nsq
├── athaliana_db.translation
```
It is not important that you understand what each file represents. However, if you are interested in verifying that all of your chromosomes were compiled, you may view the .translation file. It should look like:
```
NC_003070.9     1
NC_003071.7     2
NC_003074.8     3
NC_003075.7     4
NC_003076.8     5
NC_037304.1     6
NC_000932.1     7
```
We see that all seven chromosomes were succesffuly compiled! We are now ready to run the RepeatModeler with the script `02_repeatmodeler.sh`.
```
module load RepeatModeler/2.0.4
module load ninja/0.95 

# go to repeat directory
REPDIR=../../results/02_mask_repeats

cd ${REPDIR}

# set repdb and run repeatmodeler
REPDB=athaliana_db

RepeatModeler -threads 30 -database ${REPDB} -LTRStruct 
```
This process may run for over a day, so be patient and do not submit the job more than once! After completion of the run, there should be a directory called RM*. Let's have a look at its contents:
```
RM_150489.*/
├── consensi.fa
├── consensi.fa.classified
├── round-1
├── round-2
├── round-3
├── round-4
├── round-5
```
Per the RepeatModeler [webpage](http://www.repeatmasker.org/RepeatModeler/), we see each file as:
```
          round-1/
               sampleDB-#.fa       : The genomic sample used in this round
               sampleDB-#.fa.lfreq : The RepeatScout lmer table
               sampleDB-#.fa.rscons: The RepeatScout generated consensi
               sampleDB-#.fa.rscons.filtered : The simple repeat/low
                                               complexity filtered
                                               version of *.rscons
               consensi.fa         : The final consensi db for this round
               family-#-cons.html  : A visualization of the model
                                     refinement process.  This can be opened
                                     in web browsers that support zooming.
                                     ( such as firefox ).
                                     This is used to track down problems
                                     with the Refiner.pl
               index.html          : A HTML index to all the family-#-cons.html
                                     files.
          round-2/
               sampleDB-#.fa       : The genomic sample used in this round
               msps.out            : The output of the sample all-vs-all
                                     comparison
               summary/            : The RECON output directory
                    eles           : The RECON family output
               consensi.fa         : Same as above
               family-#-cons.html  : Same as above
               index.html          : Same as above
          round-3/
               Same as round-2
           ..
          round-n/
```
We see that we have information about the genomic sample used in each round, a consensus seqeuence frequency matrix for the genomic sample, the generated predicted consensus sequences, and visualizations. This format is repeated for various rounds with summaries of all rounds compiled in the summary directories. Our complete, predicted consensus sequences may be found in the various "consensi" fastas. Now that we have generated our consensus sequences, we are ready to mask our genome using the RepeatMasker.  


### Masking Regions of Genomic Repetition with RepeatMasker
Now that we have identified our consensus sequences, we are ready to mask them using the RepeatMasker. RepeatMasker requires two arguments, a library of repetitive regions for your organism and the genome fasta for your organism. RepeatMasker will align the repetitive regions to your genome followed by masking those repetitive regions within your genome appropriately. Let's have a look at the RepeatMasker options:
```
RepeatMasker
::small preview of options::
   -lib
        Rather than use a database, use your own RepeatModeler consensus fasta to ammend your genome
   -small
       Returns complete .masked sequence in lower case
   -xsmall
       Returns repetitive regions in lowercase (rest capitals) rather than
       masked
   -x  Returns repetitive regions masked with Xs rather than Ns
```
We want to softmask only repetitive regions, so we will be using the option "xsmall". The complete slurm script is called `03_repeatmasker.sh`:

```bash
# load software
module load RepeatMasker/4.1.2

# set variables
REPLIB=../../results/02_mask_repeats/athaliana_db-families.fa
OUTDIR=../../results/02_mask_repeats/repeatmasker
    mkdir -p ${OUTDIR}

GENOME=../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna

RepeatMasker -dir repeatmasker_out -pa 8 -lib ${REPLIB} -gff -a -noisy -xsmall ${GENOME}
```

This will produce the following files:
```
/results/02_mask_repeats/repeatmasker/repeatmasker_out
├── GCF_000001735.4_TAIR10.1_genomic.fna.align
├── GCF_000001735.4_TAIR10.1_genomic.fna.cat.gz
├── GCF_000001735.4_TAIR10.1_genomic.fna.masked
├── GCF_000001735.4_TAIR10.1_genomic.fna.out
├── GCF_000001735.4_TAIR10.1_genomic.fna.out.gff
└── GCF_000001735.4_TAIR10.1_genomic.fna.tbl
```
We are mainly interested in the masked fasta, let's give it a quick look on `GCF_000001735.4_TAIR10.1_genomic.fna.masked`, which shows the genome is soft-masked.
```
>less  GCF_000001735.4_TAIR10.1_genomic.fna.masked
>Chr1
ccctaaaccctaaaccctaaaccctaaacctctgaatccttaatccctaa
atccctaaatctttaaatcctacatccatgaatccctaaatacctaattc
cctaaacccgaaaccGGTTTCTCTGGTTGAAAATCATTGTGTATATAATG
ATAATTTTATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTT
AAAAATATCATTTGAGGTCAATACAAATCCTATTTCTTGTGGTTTTCTTT
CCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAA
GCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTATTGTAACCTTA
GGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACTGTTTATGTTT
GGACATTTATTGTCATTCTTACTCCTTTGTGGAAATGTTTGTTCTATCAA
```
RepeatMasker produces masking stats and other relevant output files, lets have a look at few of them.

**GCF_000001735.4_TAIR10.1_genomic.fna.tbl** have summary stats showing % of genome masked and the repeat elements that were used in the masking and there individual contribution in masking.
```
>head -12 GCF_000001735.4_TAIR10.1_genomic.fna.tbl
==================================================
file name: GCF_000001735.4_TAIR10.1_genomic.fna
sequences:             7
total length:  119668634 bp  (119483030 bp excl N/X-runs)
GC level:         36.06 %
bases masked:   20608830 bp ( 17.22 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements         9948      8538996 bp    7.14 %
   SINEs:              181        23261 bp    0.02 %

```
**GCF_000001735.4_TAIR10.1_genomic.fna.out** This file provide some details on each individual masking by providing SW (Smith-Waterman) scores, percent divergence, deletion, insertion, genomic location and repeat type. SW and divergence score cutoff can bet set while running RepeatMasker.
```
 SW   perc perc perc  query     position in query              matching           repeat            position in repeat
 score   div. del. ins.  sequence  begin    end          (left)   repeat             class/family  begin   end    (left)     ID
   349   13.6  5.2  4.3  Chr1             1      115 (30427556) C rnd-1_family-2     Unspecified    (1084)    286     171     1
    21    2.9  5.7  0.0  Chr1          1064     1098 (30426573) + (CACCCCC)n         Simple_repeat       1     37     (0)     2 *
    22   10.0  0.0  0.0  Chr1          1066     1097 (30426574) + (C)n               Simple_repeat       1     32     (0)     3
    15   17.1  0.0  0.0  Chr1          1155     1187 (30426484) + (TTTCTT)n          Simple_repeat       1     33     (0)     4
    28    8.4  0.0  0.0  Chr1          4291     4328 (30423343) + (AT)n              Simple_repeat       1     38     (0)     5
    16    9.3  0.0  0.0  Chr1          5680     5702 (30421969) + (T)n               Simple_repeat       1     23     (0)     6
    36    0.0  0.0  0.0  Chr1          8669     8699 (30418972) + (CT)n              Simple_repeat       1     31     (0)     7
```
***GCF_000001735.4_TAIR10.1_genomic.fna.out.gff** provides the masking information in a gff format. The columns are , chromosome, software used to annotate the feature (here RepeatMaasker), Feature type, Start and End of feature, Score, Strand and last column shows additional attributes assosiated with the feature.
```
##gff-version 2
##date 2021-07-27
##sequence-region Athaliana_167_TAIR9.fa
Chr1    RepeatMasker    similarity      1       115     13.6    -       .       Target "Motif:rnd-1_family-2" 171 286
Chr1    RepeatMasker    similarity      1064    1098     2.9    +       .       Target "Motif:(CACCCCC)n" 1 37
Chr1    RepeatMasker    similarity      1066    1097    10.0    +       .       Target "Motif:(C)n" 1 32
Chr1    RepeatMasker    similarity      1155    1187    17.1    +       .       Target "Motif:(TTTCTT)n" 1 33
Chr1    RepeatMasker    similarity      4291    4328     8.4    +       .       Target "Motif:(AT)n" 1 38
Chr1    RepeatMasker    similarity      5680    5702     9.3    +       .       Target "Motif:(T)n" 1 23
Chr1    RepeatMasker    similarity      8669    8699     0.0    +       .       Target "Motif:(CT)n" 1 31
```
**GCF_000001735.4_TAIR10.1_genomic.fna.align** shows the alignment of repeat to the genomic sequence.
```
349 13.63 5.17 4.35 Chr1 1 115 (30427556) C rnd-1_family-2#Unspecified (1084) 286 171 m_b1s001i0 1
  Chr1                   1 CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAA 50
                                                         -  i  ii     -
C rnd-1_family-        286 CCCTAAACCCTAAACCCTAAACCCTAAACC-CTAAACTCTTAA-CCCTAA 239
  Chr1                  51 ATCCCTAAATC-TTTAAATCCTACATCCATGAATCCCTAAAT-----ACC 94
                            -       i -ii    i    v i  - i  -       ?-----
C rnd-1_family-        238 A-CCCTAAACCGCCTAAACCCTAAACCC-TAAA-CCCTAAANCCTAAACC 192
  Chr1                  95 TAATTCCCTAAACCCGAAACC 115
                           i  vv          v
C rnd-1_family-        191 CAAAACCCTAAACCCTAAACC 171
Matrix = 20p35g.matrix
Kimura (with divCpGMod) = 14.34
Transitions / transversions = 2.50 (10/4)
Gap_init rate = 0.06 (7 / 114), avg. gap size = 1.57 (11 / 7)
21 2.94 5.71 0.00 Chr1 1064 1098 (30426573) (CACCCCC)n#Simple_repeat 1 37 (0) m_b1s252i0 2
  Chr1                1064 CACCCCCCACCTCCC-CCCCCC-CCCCCCACCCCCCA 1098
                                      i   -      -
  (CACCCCC)n#Si          1 CACCCCCCACCCCCCACCCCCCACCCCCCACCCCCCA 37
Matrix = Unknown
Transitions / transversions = 1.00 (1/0)
Gap_init rate = 0.06 (2 / 34), avg. gap size = 1.00 (2 / 2)
```

### Evaluating using BUSCO  
In here we will evalute the assemblies using BUSCO.  
```bash
module load busco/5.4.5

# set output directory for busco, cd into parent directory
BUSCODIR=../../results/02_mask_repeats/busco/
mkdir -p ${BUSCODIR}
cd ${BUSCODIR}

OUTDIR=masked_genome

# input masked genome
MASKED_GENOME=../../02_mask_repeats/repeatmasker/repeatmasker_out/GCF_000001735.4_TAIR10.1_genomic.fna.masked

# green plant busco database, already downloaded on xanadu. see busco documentation for how to obtain
BUSCODB=/isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10

# run busco
busco -i ${MASKED_GENOME} \
        -o ${OUTDIR} \
        -c 8 \
        -f \
        -l ${BUSCODB} \
        -m genome
```  

General useage of the command: 
```
usage: busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS] 
``` 

The command options we will be using: 
```
-i FASTA FILE   Input sequence file in FASTA format
-l LINEAGE      Specify the name of the BUSCO lineage
-o OUTPUT       Output folders and files will be labelled with this name
-m MODE         BUSCO analysis mode
					- geno or genome, for genome assemblies (DNA)
					- tran or transcriptome, for transcriptome assemblies (DNA)
					- prot or proteins, for annotated gene sets (protein)
```
The complete BUSCO scrip is called 04_busco.sh 

In here the busco metrics proposed to describe genome/gene-set/transcriptome completeness used the the following notation: 
Where recovered genes are marked as complete (C), and complete genes found with more than one copy is depected as duplicate (D), and complete single copy genes as single-copy (S) genes. The partial recovered genes are named as fragmented (F) and the genes which could not be found is named as missing (M) genes.  
So the following notation is used in busco notation:  
C: Complete, D: duplicated, F: fragmented, M: missing, n: number of genes used.   

Summary (found at `/results/02_mask_repeats/busco/masked_genome/short_summary.specific.viridiplantae_odb10.masked_genome.txt`) of the inital assembly assesment using BUSCO:  
```
        --------------------------------------------------
        |Results from dataset viridiplantae_odb10         |
        --------------------------------------------------
        |C:99.3%[S:98.6%,D:0.7%],F:0.0%,M:0.7%,n:425      |
        |422    Complete BUSCOs (C)                       |
        |419    Complete and single-copy BUSCOs (S)       |
        |3      Complete and duplicated BUSCOs (D)        |
        |0      Fragmented BUSCOs (F)                     |
        |3      Missing BUSCOs (M)                        |
        |425    Total BUSCO groups searched               |
        --------------------------------------------------
```   

## 4. Helixer  

### Annotation Method 1 (no evidence): Helixer   
In this section we will annotate our genome with [Helixer](https://github.com/weberlab-hhu/Helixer), a program that takes only a genome as input and then performs gene calling with Deep Neural Networks.

We have pulled helixer in a container on our Xanadu, which lives in a public directory. If you are working on a different cluster you will want to pull the image yourself from docker. On Xanadu we have this container in a global location: `/isg/shared/databases/nfx_singularity_cache/helixer-docker_helixer_v0.3.0a0_cuda_11.2.0-cudnn8.sif` 

Navigate the the 03_helixer directory in `/scripts`. Helixer can be run with the following script `01_helixer.sh`

```bash
#!/bin/bash
#SBATCH --job-name=helixer
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --partition=gpu
#SBATCH --qos=general
#SBATCH --constraint="AVX2&FMA3"
#SBATCH --mail-type=ALL
#SBATCH --mem=20G
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load singularity/3.9.2
mkdir tmp

export TMPDIR=$PWD/tmp
OUTDIR=../../results/03_helixer/
mkdir -p ${OUTDIR}

GENOME=../../results/02_mask_repeats/repeatmasker/repeatmasker_out/GCF_000001735.4_TAIR10.1_genomic.fna.masked

singularity exec /isg/shared/databases/nfx_singularity_cache/helixer-docker_helixer_v0.3.0a0_cuda_11.2.0-cudnn8.sif Helixer.py --lineage land_plant --fasta-path ${GENOME} --gff-output-path ${OUTDIR}/Arabidopsis_thaliana.gff3
```
In this script we execute the container with `singularity exec`, run Helixer with command `Helixer.py`, specify the taxonomic lineage with `--lineage land_plant`, the genome with `--fasta-path`, and the output path/file with `--gff-output-path`


## 5. Helixer Evaluation

Next we moved on to evaluating the annotation with AGAT, BUSCO, and EnTAP.

[AGAT](https://agat.readthedocs.io/en/latest/) is a toolkit that has scripts to check, fix, and evaluate GFF/GTF files. On Xanadu we have agat as a container in a global location: `/isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img`, however if you are working on a different cluster you will need to download the software or pull the container.

Still in the 03_helixer directory, the first evaluation script is `02_agat_stats.sh`:

```bash
module load singularity

OUTDIR="../../results/03_helixer/agat/statistics"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sp_statistics.pl \
        --gff ../../Arabidopsis_thaliana.gff3 \
        -g ../../../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna \
        -o stats.txt


singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sq_stat_basic.pl \
        --gff ../../Arabidopsis_thaliana.gff3 \
        -g ../../../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna \
        -o basic_stats.txt

singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sp_functional_statistics.pl \
        --gff ../../Arabidopsis_thaliana.gff3 \
        -g ../../../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna \
        -o functional_statistics
```
We utilize three of AGAT's evaluation scripts here: `agat_sp_statistics.pl`, `agat_sq_stat_basic.pl`, `agat_sp_functional_statistics.pl`

Results from these runs can be found in the directory `/results/03_helixer/agat/statisticsls`. 

Let's look at the `stats.txt` (`less stats.txt`)

Here we can look at the mono/multi exonic gene ratio:

"Number of single exon gene" / ("Number of gene" - "Number of single exon gene")

5067/22015=0.230


Next we extract a single transcript per gene model and then extract the coding sequences in the script `03_agat_codingsequences.sh` in order to evaluate them with BUSCO and EnTAP. This can be done with `agat_sp_keep_longest_isoform.pl` and `agat_sp_extract_sequences.pl`.


```bash
module load singularity

OUTDIR="../../results/03_helixer/agat/codingsequences"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sp_keep_longest_isoform.pl \
        --gff ../../Arabidopsis_thaliana.gff3 \ 
        -o helixer_longest_isoform.gtf 


singularity exec /isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.0.0--pl5321hdfd78af_0.img agat_sp_extract_sequences.pl \
        -g helixer_longest_isoform.gtf \
        -f ../../../../data/genome/GCF_000001735.4_TAIR10.1_genomic.fna \
        -p \
        -o helixer_proteins.fasta.faa
```
The protein file (`helixer_proteins.fasta.faa`) can then be used in further evaluations. Next we'll run the script `04_busco.sh`:

```bash
module load busco/5.4.5

# green plant busco database. already on the xanadu cluster. see documentation for how to obtain
BUSCODB=/isg/shared/databases/BUSCO/odb10/lineages/viridiplantae_odb10

#Output directory
OUTDIR=../../results/03_helixer/busco
mkdir -p ${OUTDIR}



# predicted proteins
PROTEINS=../../results/03_helixer/agat/codingsequences/helixer_proteins.fasta.faa


# run busco
busco -i ${PROTEINS} \
        -o ${OUTDIR} \
        -c 8 \
        -l ${BUSCODB} \
        -m prot \
        -f
```
Let's take a look at the output (`/results/03_helixer/busco/short_summary.specific.viridiplantae_odb10.busco.txt`):

```bash
 ***** Results: *****

        C:99.1%[S:98.4%,D:0.7%],F:0.2%,M:0.7%,n:425        
        421     Complete BUSCOs (C)                        
        418     Complete and single-copy BUSCOs (S)        
        3       Complete and duplicated BUSCOs (D)         
        1       Fragmented BUSCOs (F)                      
        3       Missing BUSCOs (M)                         
        425     Total BUSCO groups searched


```

Finally, we will run `05_EnTAP.sh` to generate a funtional annotation:

```bash
module load EnTAP/0.9.0-beta
module load diamond/0.9.36

OUTDIR=../../results/03_helixer/EnTAP
    mkdir -p ${OUTDIR}
    
    cp entap_config.txt ${OUTDIR}
    cd ${OUTDIR}

PROTEINS=../agat/codingsequences/helixer_proteins.fasta.faa

EnTAP --runP \
        -i ${PROTEINS} \
        -d /isg/shared/databases/Diamond/RefSeq/complete.protein.faa.216.dmnd \
        -d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd \
        --ontology 0 \
        --threads 16 

date
```
The output will be writen in the the entap_outfiles/ folder.

```bash
entap_outfiles
├── final_results
├── ontology
├── similarity_search
└── transcriptomes
```

Similarity search directory will contain the results from the [diamond](https://github.com/bbuchfink/diamond) database(s) you included in your search. Inside the processed/ folder you will find the information based on the hits returned from similarity searching against each database. This information contains the best hits (discussed previously) from each database based on e-value, coverage, informativeness, phylogenetic closeness, and contaminant status.

- In the database you selected under the processed directory it will contain best_hits.faa and .fnn and .tsv files:
- best_hits_contam.faa/.fnn/.tsv will contain contaminants (protein/nucleotide) separated from the best hits file.
- best_hits_no_contam.faa/.fnn/.tsv will contain sequences (protein/nucleotide) that were selected as best hits and not flagged as contaminants
- no_hits.faa/.fnn/.tsv contain sequences (protein/nucleotide) from the transcriptome that did not hit against this particular database
- unselected.tsv will contain result in several hits for each query sequence. With only one best hit being selected, the rest are unselected and end up here

In the Ontology search folder we will find the ortholog groups / ontolgoy results agains EggNOG pipeline. Inside the processed folder it will contain the annotations.

In the final_results folder, it will contain the final EnTAP annotations. These files are the summation of each stage of the pipeline and contain the combined information. So these can be considered the most important files! Gene ontology terms are normalized to levels based on the input flag from the user (or the default of 0,3,4). A level of 0 within the filename indicates that ALL GO terms will be printed to the annotation file. Normalization of GO terms to levels is generally done before enrichment analysis and is based upon the hierarchical setup of the Gene Ontology database.

- final_annotations_lvlX.tsv : X represents the normalized GO terms for the annotation
- final_annotated.faa / .fnn : Nucleotide and protein fasta files containing all sequences that either hit databases through similarity searching or through the ontology stage
- final_unannotated.aa / .fnn : Nucleotide and protein fasta files containing all sequences that did not hit either through similarity searching nor through the ontology stage

In the following example it shows the results agains the complete diamond database we selected, and EggNOG pipeline and the final results.

More information on EnTAP can be found in [EnTAP](https://entap.readthedocs.io/) documentation, which has a very comprehensive description.

Let's take a look at some results (`/results/03_helixer/EnTAP/entap_outfiles/log_file_2023Y6M6D-13h14m4s.txt`)

```bash
Final Annotation Statistics
------------------------------------------------------
Total Sequences: 27082
Similarity Search
        Total unique sequences with an alignment: 25555
        Total unique sequences without an alignment: 1527
Gene Families
        Total unique sequences with family assignment: 25834
        Total unique sequences without family assignment: 1248
        Total unique sequences with at least one GO term: 20969
        Total unique sequences with at least one pathway (KEGG) assignment: 6038
Totals
        Total unique sequences annotated (similarity search alignments only): 411
        Total unique sequences annotated (gene family assignment only): 690
        Total unique sequences annotated (gene family and/or similarity search): 26245
        Total unique sequences unannotated (gene family and/or similarity search): 837
```

## 6. Braker Annotation with Protein Evidence

Our second method of annotation in this tutorial is using protein evidence as input with Braker. First we must obtain the protein data and concatenate it into a single fasta file. Navigate to the `scripts/04_braker` directory, and run the `/01_get_proteins.sh` script:

```bash
OUTDIR=../../results/04_braker/proteins
    mkdir -p ${OUTDIR}
    cd ${OUTDIR}

# per the braker documentation, a protein database can be obtained as below
    # this code downloads protein sequences in fasta format for all genes in orthoDB for the Viridiplantae
    # then it concatenates them into a single fasta file
    # https://github.com/gatech-genemark/ProtHint#protein-database-preparation
wget --no-check-certificate https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz
tar -xvzf odb10_plants_fasta.tar.gz
cat plants/Rawdata/*fs >proteins.fa
```
Now we can run Braker with the proteins and masked genome as input (`02_braker_proteins.sh`):

```bash
# load software
module load BRAKER/2.1.6
module load ProtHint/2.6.0

# for testing new perl version:
module unload perl
module load perl/5.36.0

# set braker output directory and cd there
OUTDIR=../../results/04_braker/proteins
    mkdir -p ${OUTDIR}
cd ${OUTDIR}

# copy augustus config, set variable
export AUGUSTUS_CONFIG_PATH=$(pwd)/config
    cp -r /isg/shared/apps/augustus/3.6.0/config/ config

# create a temp directory for braker temp files
export TMPDIR=$(pwd)/braker_tmp 
    mkdir -p ${TMPDIR}

# set variables for input/output files
GENOME="../../02_mask_repeats/repeatmasker/repeatmasker_out/GCF_000001735.4_TAIR10.1_genomic.fna.masked"
PROTEINDB=proteins.fa

# run braker
braker.pl \
    --genome=${GENOME} \
    --prot_seq=${PROTEINDB} \
    --softmasking \
    --cores 16 \
    --gff3
```
To evaluate this genome we will follow the same steps we did with Helixer (agat, busco, and EnTAP). Let's have a look at some of the results:

AGAT:

Mono/Multi Ratio
"Number of single exon gene" / ("Number of gene" - "Number of single exon gene")
7132/23114=0.309

Busco:
```bash
	***** Results: *****

	C:99.5%[S:98.8%,D:0.7%],F:0.2%,M:0.3%,n:425	   
	423	Complete BUSCOs (C)			   
	420	Complete and single-copy BUSCOs (S)	   
	3	Complete and duplicated BUSCOs (D)	   
	1	Fragmented BUSCOs (F)			   
	1	Missing BUSCOs (M)			   
	425	Total BUSCO groups searched		   

```

EnTAP:
 ```bash
 Final Annotation Statistics
------------------------------------------------------
Total Sequences: 30246
Similarity Search
        Total unique sequences with an alignment: 27785
        Total unique sequences without an alignment: 2461
Gene Families
        Total unique sequences with family assignment: 28008
        Total unique sequences without family assignment: 2238
        Total unique sequences with at least one GO term: 22365
        Total unique sequences with at least one pathway (KEGG) assignment: 6223
Totals
        Total unique sequences annotated (similarity search alignments only): 925
        Total unique sequences annotated (gene family assignment only): 1148
        Total unique sequences annotated (gene family and/or similarity search): 28933
        Total unique sequences unannotated (gene family and/or similarity search): 1313
        
 ```

Method 3: EASEL

EASEL (Efficient, Accurate, Scalable Eukaryotic modeLs) is a genome annotation pipeline that utilizes machine learning, RNA folding, and functional annotations to enhance gene prediction accuracy. is a nextflow pipeline that annotates a genome with RNA-seq reads as evidence. 

More about methods and the workflows may be found on the EASEL website[https://gitlab.com/PlantGenomicsLab/easel].

To run EASEL on xanadu we must first pull the container, and then provide it with a set of paramaters. The params.yaml file:

Then we ran the following script:


## 19. GFFCompare   
GFFCompare can be used to compare, merge, annotate and estimate accuracy of one or more GFF files when compared with a reference annotation. Here we run GFFCompare with the three annotation files we generated through this tutorial, and the "truth" that we pulled in the beginning in the 01_raw_data directory. 

Navigate to `06_gffcompare` and `sbatch 01_gffcompare.sh`:

```bash
module load gffcompare/0.10.4

# output directory
OUTDIR=../../results/06_gffcompare_2
    mkdir -p ${OUTDIR}

# files
NCBIANNO=../../data/genome/GCF_000001735.4_TAIR10.1_genomic.gff
HELIXER=../../results/03_helixer/Arabidopsis_thaliana.gff3
BRAKERPROTEIN=../../results/04_braker/proteins/braker/augustus.hints.gtf
EASEL=../../results/05_EASEL/arabidopsis/final_predictions/arabidopsis_filtered.gtf

# run gffcompare
gffcompare -R -r <(awk '$7 !~ /?/' ${NCBIANNO}) -o ${OUTDIR}/helixer ${HELIXER}

gffcompare -R -r <(awk '$7 !~ /?/' ${NCBIANNO}) -o ${OUTDIR}/braker ${BRAKERPROTEIN}

gffcompare -R -r <(awk '$7 !~ /?/' ${NCBIANNO}) -o ${OUTDIR}/easel ${EASEL}
```

Let's take a look at the results: 

EASEL:
```bash
#= Summary for dataset: ../../results/05_EASEL/arabidopsis/final_predictions/arabidopsis_filtered.gtf 
#     Query mRNAs :   29381 in   19324 loci  (25053 multi-exon transcripts)
#            (6209 multi-transcript loci, ~1.5 transcripts per locus)
# Reference mRNAs :   38462 in   20342 loci  (34229 multi-exon)
# Super-loci w/ reference transcripts:    18924
#-----------------| Sensitivity | Precision  |
        Base level:    68.0     |    99.2    |
        Exon level:    52.9     |    65.9    |
      Intron level:    77.9     |    83.3    |
Intron chain level:    29.9     |    40.8    |
  Transcript level:    29.6     |    38.7    |
       Locus level:    54.6     |    57.1    |

     Matching intron chains:   10225
       Matching transcripts:   11366
              Matching loci:   11112

          Missed exons:   18649/159922	( 11.7%)
           Novel exons:     942/126977	(  0.7%)
        Missed introns:   11149/115825	(  9.6%)
         Novel introns:    1786/108273	(  1.6%)
           Missed loci:       0/20342	(  0.0%)
            Novel loci:      15/19324	(  0.1%)

 Total union super-loci across all input datasets: 18946 
29381 out of 29381 consensus transcripts written in ../../results/06_gffcompare_2/easel.annotated.gtf (0 discarded as redundant)
```
BRAKER:
```bash
#= Summary for dataset: ../../results/04_braker/proteins/braker/augustus.hints.gtf 
#     Query mRNAs :   32230 in   30204 loci  (24877 multi-exon transcripts)
#            (1655 multi-transcript loci, ~1.1 transcripts per locus)
# Reference mRNAs :   48521 in   27544 loci  (41185 multi-exon)
# Super-loci w/ reference transcripts:    27253
#-----------------| Sensitivity | Precision  |
        Base level:    68.8     |    98.0    |
        Exon level:    56.1     |    68.0    |
      Intron level:    84.7     |    92.3    |
Intron chain level:    39.8     |    65.9    |
  Transcript level:    37.8     |    56.9    |
       Locus level:    64.8     |    59.7    |

     Matching intron chains:   16405
       Matching transcripts:   18327
              Matching loci:   17838

          Missed exons:   15214/186166	(  8.2%)
           Novel exons:    3582/151634	(  2.4%)
        Missed introns:   12310/131340	(  9.4%)
         Novel introns:    6505/120548	(  5.4%)
           Missed loci:       0/27544	(  0.0%)
            Novel loci:     926/30204	(  3.1%)

 Total union super-loci across all input datasets: 28276 
32230 out of 32230 consensus transcripts written in ../../results/06_gffcompare_2/braker.annotated.gtf (0 discarded as redundant)
```
HELIXER:
```bash
#= Summary for dataset: ../../results/03_helixer/Arabidopsis_thaliana.gff3 
#     Query mRNAs :   27082 in   27082 loci  (22015 multi-exon transcripts)
#            (0 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   46994 in   26068 loci  (40657 multi-exon)
# Super-loci w/ reference transcripts:    24454
#-----------------| Sensitivity | Precision  |
        Base level:    86.7     |    95.7    |
        Exon level:    70.4     |    80.3    |
      Intron level:    84.5     |    89.0    |
Intron chain level:    34.4     |    63.6    |
  Transcript level:    36.5     |    63.3    |
       Locus level:    65.1     |    63.3    |

     Matching intron chains:   13999
       Matching transcripts:   17147
              Matching loci:   16981

          Missed exons:    5266/183080	(  2.9%)
           Novel exons:    5195/150243	(  3.5%)
        Missed introns:    5319/129769	(  4.1%)
         Novel introns:    5691/123161	(  4.6%)
           Missed loci:       0/26068	(  0.0%)
            Novel loci:     514/27082	(  1.9%)

 Total union super-loci across all input datasets: 25087 
27082 out of 27082 consensus transcripts written in ../../results/06_gffcompare_2/helixer.annotated.gtf (0 discarded as redundant)
```

 
