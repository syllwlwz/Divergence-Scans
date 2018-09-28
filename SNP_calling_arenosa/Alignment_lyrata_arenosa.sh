#!/bin/bash

#Alignment
#samtools 1.3.1 with htslib 1.3.1
#bwa 0.7.15

#replace number with line number in Cutadapt_arenosa.list
file=$(cat Cutadapt_arenosa.list | head -n number| tail -n 1)

./software/bwa-0.7.15/bwa mem -t 12 -k 10 alygenomes.fasta -M cutadapt_arenosa/$file.R1.cutadapt.fastq.gz cutadapt_arenosa/$file.R2.cutadapt.fastq.gz > aligned_arenosa/$file.sam 2> aligned_arenosa/$file.sam.err

#t: number of threads
#k: minimum length of seed region
 
./software/samtools-1.3/samtools view -b aligned_arenosa/$file.sam | ./software/samtools-1.3/samtools sort -T $file > aligned_arenosa/$file.sort.bam 2> aligned_arenosa/$file.sort.bam.err

#b: output in bam format
#T: write temporary files

./software/samtools-1.3/samtools index aligned_arenosa/$file.sort.bam

./software/samtools-1.3/samtools idxstats aligned_arenosa/$file.sort.bam > aligned_arenosa/$file.idxstats

./software/samtools-1.3/samtools flagstat aligned_arenosa/$file.sort.bam > aligned_arenosa/$file.flagstat



