#!/bin/bash

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l select=1:ncpus=4:mem=15gb
#PBS -l pvmem=12gb
#PBS -l place=pack:shared
#PBS -l walltime=24:00:00
#PBS -l cput=24:00:00
#PBS -M jamesthornton@email.arizona.edu
#PBS -m bea

while read FASTA; do
    FILENAME=$(basename $FASTA | cut -d '.' -f 1)
    $DIAMOND_DIR/diamond blastx -p $THREADS -d $DIAMOND_DB -q $FASTA -o $OUT_DIR/$FILENAME.dout    
done < $FASTA_LIST
