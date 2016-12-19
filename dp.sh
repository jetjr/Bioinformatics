#!/bin/bash

DB=../neutropenicfever/swissprot2taxonomy
IDS=../neutropenicfever/diamond_out/IonXpress_001_NF002.dout.ids
OUT=out.txt

cat /dev/null > $OUT

i=0
while read ID; do
  let i++
  printf "%3d: %s\n" $i $ID
  awk -F"\t" \'"\$2 == $ID { print \$3 }"\' $DB
  #awk -F"\t" "\$2 == $ID { print \$3 }" $DB >> $OUT
  break
done < $IDS

echo Done.
