#!/bin/bash

export IDS=$1
export MAP=$2

for id in $(cat $IDS); do
    grep "$id" $MAP | awk ' BEGIN { FS="\t" } /1/ { print $2 }'
done >> pathways.txt

cat pathways.txt | sort | uniq -c > pathway_counts.txt
