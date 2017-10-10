#!/bin/bash

#PBS -W group_list=bhurwitz
#PBS -q windfall
#PBS -l jobtype=cluster_only
#PBS -l select=1:ncpus=12:mem=23gb
#PBS -l pvmem=22gb
#PBS -l walltime=24:00:00
#PBS -l cput=24:00:00
#PBS -M jamesthornton@email.arizona.edu
#PBS -m bea

#--------------EDIT THESE---------------
FASTA_DIR="/rsgrps/bhurwitz/jetjr/neutropenicfever/fasta"
OUT_DIR="/rsgrps/bhurwitz/jetjr/neutropenicfever/centrifuge"
#---------------------------------------

CENT_DB="/rsgrps/bh_class/b_compressed+h+v/b_compressed+h+v"

cd "$FASTA_DIR"
export FASTA_LIST="$FASTA_DIR/fasta-list"
ls *.fa > $FASTA_LIST
echo "FASTA files to be processed:" $(cat $FASTA_LIST)

while read FASTA; do

FILE_NAME=$(basename $FASTA | cut -d '.' -f 1)

centrifuge -x $CENT_DB -U $FASTA -S $OUT_DIR/$FILE_NAME-centrifuge_hits.tsv --report-file $OUT_DIR/$FILE_NAME-centrifuge_report.tsv -f

done < $FASTA_LIST
