#!/bin/bash

#CENTRIFUGE HIT FILE FULL PATH
export CENT_HIT="/rsgrps/bhurwitz/jetjr/neutropenicfever/final/centrifuge/NF003-centrifuge_hits.tsv"

#FASTA FILE TO EXTRACT SEQS
export FASTA="/rsgrps/bhurwitz/jetjr/neutropenicfever/final/unmapped/qced/NF003.unmapped"

#TAXONOMY ID TO EXTRACT
export ID="277"

#OUT DIRECTORY TO STORE SEQ IDS AND SEQUENCES
export OUT="/rsgrps/bhurwitz/jetjr"

./centrifuge_tax_parse.R -f $CENT_HIT -i $ID -o $OUT

module load samtools/1/1.3.1

FASTA_BASE=$(basename $FASTA | cut -d '.' -f 1)
HIT_BASE=$(basename $CENT_HIT | cut -d '-' -f 1)

#echo $FILE_BASE
xargs samtools faidx $FASTA < $OUT/$HIT_BASE-$ID.ids > $OUT/$FASTA_BASE-$ID.seqs
#
#done
