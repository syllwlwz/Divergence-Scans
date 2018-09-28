#!/bin/bash

#Namefix
#picard 2.7.1

#replace number with line number in Cutadapt_arenosa.list
file=$(cat Cutadapt_arenosa.list | head -n number | tail -n 1)

java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=dedup_arenosa/$file.dedup.bam O=namefixed_arenosa/$file.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=$file RGPL=illumina RGPU=AAAAAA RGSM=$file 2> namefixed_arenosa/$file.err
