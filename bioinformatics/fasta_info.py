#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq

total_len = 0
total_seqs = 0
len_list = []
A = 0
T = 0
C = 0
G = 0
N = 0

for seq_record in SeqIO.parse("virus_contigs.fasta", "fasta"):
    total_seqs = total_seqs + 1
    total_len = len(str(seq_record.seq)) + total_len
    len_list.append(len(str(seq_record.seq)))
    A = Seq(str(seq_record.seq)).count("A") + A
    T = Seq(str(seq_record.seq)).count("T") + T
    C = Seq(str(seq_record.seq)).count("C") + C
    G = Seq(str(seq_record.seq)).count("G") + G
    N = Seq(str(seq_record.seq)).count("N") + N

large = str(max(len_list))
small = str(min(len_list))
    
print("Total Sequences: ","{0:>4}".format(total_seqs))
print("Largest Sequence:","{0:>7}".format(large))
print("Smallest Sequence: ","{0:>3}".format(small))

print("-"*25)

print("Adenines: ","{0:>5}".format(A))
print("Thymines: ","{0:>5}".format(T))
print("Cytosines:",C)
print("Guanines: ","{0:>5}".format(G))
print("Ns: ","{0:>7}".format(N))
print("Total: ","{0:>9}".format(total_len))
