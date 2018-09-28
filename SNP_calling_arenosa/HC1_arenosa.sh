#!/bin/bash

#GATK version 3.7
#cohorts separately

#replace number with line number in Cutadapt_arenosa.list
file=$(cat Cutadapt_arenosa.list | head -n $SGE_TASK_ID | tail -n 1)

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T HaplotypeCaller -I realigned_arenosa/$file.realigned.bam --min_base_quality_score 25 --min_mapping_quality_score 25 -rf DuplicateRead -rf BadMate -rf BadCigar -ERC BP_Resolution -R alygenomes.fasta  -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L scaffold_6 -L scaffold_7 -L scaffold_8 -o HC1_arenosa/$file.vcf.gz -ploidy 4 --pcr_indel_model NONE -nct 4 --heterozygosity 0.04 -rf NotPrimaryAlignment --indel_heterozygosity 0.02 2> HC1_arenosa/$file.err
