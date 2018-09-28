#!/bin/bash

#Realignment
#GATK3.7

#run from pflaphy-gscan

#replace number with line number in Cutadapt_arenosa.list
file=$(cat Cutadapt_arenosa.list | head -n number| tail -n 1)

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T RealignerTargetCreator -R alygenomes.fasta -o namefixed_arenosa/$file.IndelRealigner.intervals -I namefixed_arenosa/$file.dedup.bam 2> realigned_arenosa/$file.Interval.err; java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T IndelRealigner -targetIntervals namefixed_arenosa/$file.IndelRealigner.intervals -I namefixed_arenosa/$file.dedup.bam -R alygenomes.fasta -o realigned_arenosa/$file.realigned.bam 2> realigned_arenosa/$file.err

