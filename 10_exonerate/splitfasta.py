#!/usr/bin/python

from Bio import SeqIO
import argparse
#import numpy as np
import os

parser = argparse.ArgumentParser(
     prog='splitfasta.py',
     usage='''python splitfasta.py --fasta [Genome fasta file] --path [Path of genome file] --pieces [No. of pieces desired] --pathOut [direct the output]''',
     description='''This program splits a fasta sequence into several similarly-sized pieces.''',
     epilog='''It requires numpy and biopython libraries''')
parser.add_argument('--fasta', type=str, help='The name of the fasta file', required=True)
parser.add_argument('--path', type=str, help='The path of the fasta file', required=False)
parser.add_argument('--pieces', type=int, help='No. of pieces desired', required=True)
parser.add_argument('--pathOut', type=str, help='path of output files', required=False)

args = parser.parse_args()
pathOut = args.pathOut
filename1=args.fasta
path=args.path
pieces=args.pieces
if path==None:
  filename2=open(filename1,'r')
else:
  filename2=open(os.path.join(path, filename1), "r")

seq_records1 = SeqIO.parse(filename2, "fasta")
seq_records1 = list(seq_records1)

# Calulate total number of bases

genome=0

for seq_record in seq_records1:
   genome+=len(seq_record.seq)

splice=genome/pieces  
j=0 
for i in range(pieces):
    if pathOut==None:
      c = filename1 + str(i+1) + ".fa"
    else:
      c = pathOut + filename1 + str(i+1) + ".fa"
    breaks=0
    with open(c,'w') as ijk:        
         while (breaks < splice) and j < len(seq_records1):
             ijk.write("%s%s\n%s\n" % (">",seq_records1[j].id, seq_records1[j].seq)) # Making the fasta file
             breaks+=len(seq_records1[j].seq)
             j+=1
    ijk.close()

