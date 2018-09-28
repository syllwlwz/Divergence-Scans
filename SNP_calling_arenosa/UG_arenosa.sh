#!/bin/bash

#GATK version 3.7

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T UnifiedGenotyper -I KletKowaarenosa_UG.list --min_base_quality_score 25 -rf DuplicateRead -rf BadMate -rf BadCigar -R alygenomes.fasta  -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L scaffold_6 -L scaffold_7 -L scaffold_8 -o UG_arenosa/KletKowa_UG.vcf.gz -ploidy 4 -nct 12 -rf NotPrimaryAlignment -glm SNP --heterozygosity 0.04 -stand_call_conf 25 --output_mode EMIT_VARIANTS_ONLY -dcov 200 2> UG_arenosa/KletKowa_UG.err

