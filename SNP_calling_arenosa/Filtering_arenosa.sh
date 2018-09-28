#!/bin/bash

#GATK3.7
#cohorts separately

#Filter 1: select bialleleic SNPs without missing genotypes
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V HC2_arenosa/MiasZapa_arenosa.vcf -R alygenomes.fasta -selectType SNP -select "AN == 60" -restrictAllelesTo BIALLELIC -env -o Filtered_arenosa/MiasZapa_arenosaFil_AN.vcf 2> Filtered_arenosa/MiasZapaFil_AN.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa/MiasZapa_arenosaFil_AN.vcf -R alygenomes.fasta -F CHROM -F POS -F SOR -F QD -F FS -F ReadPosRankSum -F MQ -F MQRankSum -raw --allowMissingData -o Filtered_arenosa/MiasZapa_Filtest.table 2> Filtered_arenosa/MiasZapa_Filtest.err

#Filter 1: select bialleleic SNPs without missing genotypes
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V HC2_arenosa/KletKowa_arenosa.vcf -R alygenomes.fasta -selectType SNP -select "AN == 68" -restrictAllelesTo BIALLELIC -env -o Filtered_arenosa/KletKowa_arenosaFil_AN.vcf 2> Filtered_arenosa/KletKowaFil_AN.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa/KletKowa_arenosaFil_AN.vcf -R alygenomes.fasta -F CHROM -F POS -F SOR -F QD -F FS -F ReadPosRankSum -F MQ -F MQRankSum -raw --allowMissingData -o Filtered_arenosa/KletKowa_Filtest.table 2> Filtered_arenosa/KletKowa_Filtest.err

