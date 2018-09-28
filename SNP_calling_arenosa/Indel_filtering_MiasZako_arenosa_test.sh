#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_MiasZako_indel
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

#ls HC1_arenosa/*.vcf.gz | grep "Mias" | grep "vcf.gz" | sed 's/.vcf.gz//g' | sed 's/HC1_arenosa\///g' > Mias.list
#ls HC1_arenosa/*.vcf.gz | grep "Zako" | grep "vcf.gz" | sed 's/.vcf.gz//g' | sed 's/HC1_arenosa\///g' > Zako.list

#Filter 1: select bialleleic SNPs without missing genotypes
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V HC2_arenosa/MiasZapa_arenosa.vcf -R alygenomes.fasta -selectType INDEL -o Filtered_arenosa/MiasZapa_INDEL_arenosa.vcf 2> Filtered_arenosa/MiasZapa_INDEL_arenosa.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa/MiasZapa_INDEL_arenosa.vcf -R alygenomes.fasta -F CHROM -F POS -F SOR -F QD -F FS -F ReadPosRankSum -raw --allowMissingData -o Filtered_arenosa/MiasZapa_INDEL_arenosa.table 2> Filtered_arenosa/MiasZapa_INDEL_arenosa2.err

#Filter 1: select bialleleic SNPs without missing genotypes
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V HC2_arenosa/KletKowa_arenosa.vcf -R alygenomes.fasta -selectType INDEL -o Filtered_arenosa/KletKowa_INDEL_arenosa.vcf 2> Filtered_arenosa/KletKowa_INDEL_arenosa.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa/KletKowa_INDEL_arenosa.vcf -R alygenomes.fasta -F CHROM -F POS -F SOR -F QD -F FS -F ReadPosRankSum -raw --allowMissingData -o Filtered_arenosa/KletKowa_INDEL_arenosa.table 2> Filtered_arenosa/KletKowa_INDEL_arenosa2.err

#Filter 1: select bialleleic SNPs without missing genotypes
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V HC2_arenosa/WulmCho_arenosa.vcf -R alygenomes.fasta -selectType INDEL -o Filtered_arenosa/WulmCho_INDEL_arenosa.vcf 2> Filtered_arenosa/WulmCho_INDEL_arenosa.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa/WulmCho_INDEL_arenosa.vcf -R alygenomes.fasta -F CHROM -F POS -F SOR -F QD -F FS -F ReadPosRankSum -raw --allowMissingData -o Filtered_arenosa/WulmCho_INDEL_arenosa.table 2> Filtered_arenosa/WulmCho_INDEL_arenosa2.err


