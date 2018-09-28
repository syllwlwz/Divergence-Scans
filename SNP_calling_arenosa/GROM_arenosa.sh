#!/bin/bash

#GROM
#samtools 1.3.1 with htslib 1.3.1
#replace number with line number in GROM_arenosa.list

file=$(cat GROM_arenosa.list | head -n number | tail -n 1)
./software/GROM/dist/GROM -i GROM_analysis/$file.bam -r alygenomes.fasta -o GROM_analysis/$file.predictions.vcf -P 6 -b 25 -q 25 -v 0.04 -e 0.02 -V 0.04 -U 3 -A 4 -p 4 -s 30 2> GROM_analysis/$file.err

