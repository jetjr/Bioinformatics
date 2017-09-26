#!/usr/bin/env python3

import math
import sys

from Bio import SeqIO

args = sys.argv

if len(args) < 3:
    print('Usage:', args[0], 'FASTA', 'Size')
    sys.exit(1)

header = 1
size = int(args[2])

for seq_record in SeqIO.parse(args[1], "fasta"):
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
  
