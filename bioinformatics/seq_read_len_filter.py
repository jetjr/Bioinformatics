#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os

args = sys.argv

if len(args) < 2:
  print("Usage:", args[0], "FASTA", "Minimum Length")
  sys.exit(1)

seq_len = 0
seq_num = 0

basename = os.path.basename(args[1])
out_file = basename.split('.')[0] + "_length_distribution.txt"

for seq_record in SeqIO.parse(args[1], "fasta"):
    seq_num = seq_num + 1
    seq_len = len(str(seq_record.seq))
    if int(seq_len) > int(args[2]):
        print(">",seq_record.description, sep = "")
        print(seq_record.seq)
