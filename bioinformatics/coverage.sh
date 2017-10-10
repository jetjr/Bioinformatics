#!/bin/bash

SEQ="/rsgrps/bhurwitz/chunan/fourthtest/assembly/s_aureus_ref.fasta"
SEQ_N=$(basename $SEQ | cut -d '.' -f 1)

SAMPLE="/rsgrps/bhurwitz/chunan/fourthtest/fastqdata/IonXpress_035_21V.fasta"

module load samtools
module load bowtie2

bowtie2-build $SEQ $SEQ_N

bowtie2 --threads 2 -x $SEQ_N -U $SAMPLE -S $SEQ_N.sam -f --very-sensitive-local

samtools view -F 4 -bS $SEQ_N.sam > $SEQ_N.bam
samtools sort -T $SEQ_N.sorted -o $SEQ_N.sorted.bam $SEQ_N.bam
samtools depth $SEQ_N.sorted.bam > $SEQ_N.coverage
