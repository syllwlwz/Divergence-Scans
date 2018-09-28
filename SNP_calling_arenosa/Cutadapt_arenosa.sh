#!/bin/bash

#replace number with line number in Cutadapt_arenosa.list
file=$(cat Cutadapt_arenosa.list | head -n number| tail -n 1)
./software/cutadapt -e 0.15 -O 4 -m 25 -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA' -o cutadapt_arenosa/$file.R1.cutadapt.fastq.gz -p cutadapt_arenosa/$file.R2.cutadapt.fastq.gz Arenosa_fastq/$file.R1.fastq.gz Arenosa_fastq/$file.R2.fastq.gz > cutadapt_arenosa/$file.cutadapt 2> cutadapt_arenosa/$file.cutadapt.err


