#!/usr/bin/env python3

import math

from Bio import SeqIO

header = 1
size = 150

for seq_record in SeqIO.parse("virus_contigs.fasta", "fasta"):
    seq_len = len(seq_record.seq)
    num = range(math.ceil(seq_len / size))
    print(seq_len)
    begin = 0
    end = 150
    for i in num:
        print(">",header)
        print(str(seq_record.seq[begin:end]))
        begin = begin + size
        end = end + size
        header = header + 1
  
