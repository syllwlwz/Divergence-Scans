#!/bin/bash

#bedtools version 2.26
#cohorts separately

#run from pflaphy-gscan

#$ -t 6
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Coverage_exons
#$ -l vf=16G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Coverage_arenosa.list | head -n $SGE_TASK_ID | tail -n 1)
./software/bedtools2-2.27.1/bin/bedtools coverage -hist -sorted -g alygenomes.fasta.fai -a Alyrata_exons_sorted.gff -b realigned_arenosa/$file*realigned.bam > Coverage/Hist_arenosa_exons_$file.bed 2> Coverage/Hist_arenosa_exons_$file.err
