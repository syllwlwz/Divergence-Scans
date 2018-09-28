#!/bin/bash

#Deduplication
#picard tools 2.7.1
#replace number with line number in Cutadapt_arenosa.list

file=$(cat Cutadapt_arenosa.list | head -n number | tail -n 1)

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/picard.jar  MarkDuplicates I=aligned_arenosa/$file.sort.bam O=dedup_arenosa/$file.dedup.bam M=dedup_arenosa/duplicateMetricsFile_$file VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> dedup_arenosa/$file.err
