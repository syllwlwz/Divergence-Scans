#!/bin/bash

#bedtools version 2.26
#replace number with line number in Coverage_arenosa.list

file=$(cat Coverage_arenosa.list | head -n number | tail -n 1)
./software/bedtools2-2.27.1/bin/bedtools coverage -hist -sorted -g alygenomes.fasta.fai -a Alyrata_genes_sorted.gff -b realigned_arenosa/$file*realigned.bam > Coverage/Hist_arenosa_$file.bed 2> Coverage/Hist_arenosa_$file.err
