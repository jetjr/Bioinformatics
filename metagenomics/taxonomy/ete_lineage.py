#!/usr/bin/env python

from ete3 import NCBITaxa
import sys
import os


args = sys.argv

if len(args) < 2:
  print("Usage:", args[0], "[IDs]")
  sys.exit(1)

ncbi = NCBITaxa()

for id in open(args[1]):
    print ncbi.get_lineage(id)
