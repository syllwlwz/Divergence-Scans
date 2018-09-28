#!/bin/bash

#GATK version 3.7
#cohorts separately

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T GenotypeGVCFs -R alygenomes.fasta -nt 12 -o HC2_arenosa/KletKowa_arenosa.vcf -hets 0.04 -newQual -V KletKowaarenosa.list --indel_heterozygosity 0.02 2> HC2_arenosa/KletKowaarenosa.err

