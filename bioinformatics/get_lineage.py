#!/usr/bin/env python3

from Bio import Entrez
import sys
import os

Entrez.email = "jamesthornton@email.arizona.edu"

args = sys.argv

if len(args) < 2:
  print("Usage:", args[0], "[IDs]")
  sys.exit(1)

for id in open(args[1]):
    handle = Entrez.efetch(db="taxonomy", id=id, retmode="xml")
    records = Entrez.read(handle)
    print(records[0]["Lineage"])
