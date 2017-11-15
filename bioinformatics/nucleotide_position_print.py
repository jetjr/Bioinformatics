#!/usr/bin/env python3

import sys

args = sys.argv

seq_file = open(args[1], 'r')

print(seq_file[0:100])
