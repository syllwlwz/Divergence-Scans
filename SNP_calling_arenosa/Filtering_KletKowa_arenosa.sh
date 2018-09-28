#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_KletKowa
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

#ls HC1_arenosa/*.vcf.gz | grep "Klet" | grep "vcf.gz" | sed 's/.vcf.gz//g' | sed 's/HC1_lyrata\///g' > Klet.list
#ls HC1_arenosa/*.vcf.gz | grep "Kowa" | grep "vcf.gz" | sed 's/.vcf.gz//g' | sed 's/HC1_lyrata\///g' > Kowa.list

#Filter 1: select bialleleic SNPs without missing genotypes
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V HC2_arenosa/KletKowa_arenosa.vcf -R alygenomes.fasta -selectType SNP -select "AN == 68" -restrictAllelesTo BIALLELIC -env -o Filtered_arenosa/KletKowaFil_AN_arenosa.vcf 2> Filtered_arenosa/KletKowaFil_AN_arenosa.err

#Filter 2: GATK best practices
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_arenosa/KletKowaFil_AN_arenosa.vcf -R alygenomes.fasta -filterName "BP" -filter "QD<2.0 || FS>30.0 || MQ<50.0 || MQRankSum< -2.5 || ReadPosRankSum < -4.0 || SOR>4.0" -o Filtered_arenosa/KletKowaFil_AN_BP_arenosa.vcf 2> Filtered_arenosa/KletKowaFil_AN_BP_arenosa.err

#Filter 3: Minimum depth of 2 in every sample
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_arenosa/KletKowaFil_AN_BP_arenosa.vcf -R alygenomes.fasta -G_filter "DP<4.0" -G_filterName "lowDP" -o Filtered_arenosa/KletKowaFil_AN_DP_arenosa.list 2> Filtered_arenosa/KletKowaFil_AN_DP_arenosa.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa/KletKowaFil_AN_BP_arenosa.vcf -R alygenomes.fasta -XL Filtered_arenosa/KletKowaFil_AN_DP_arenosa.list -o Filtered_arenosa/KletKowaFil_AN_DP_arenosa.vcf 2> Filtered_arenosa/KletKowaFil_AN_DP2_arenosa.err

grep -v "#" Filtered_arenosa/KletKowaFil_AN_DP_arenosa.list | grep "lowDP" >> Filtered_arenosa/KletKowaFil_AN_DP_arenosa.vcf 2> Filtered_arenosa/KletKowaFil_lowDP_arenosa.err

#Mask filter 3 to filter 2
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_arenosa/KletKowaFil_AN_BP_arenosa.vcf -R alygenomes.fasta --mask Filtered_arenosa/KletKowaFil_AN_DP_arenosa.vcf --maskName "DP" -o Filtered_arenosa/KletKowaFil_AN_BP_DP_arenosa.vcf 2> Filtered_arenosa/KletKowaFil_AN_BP_DP_arenosa.err

#Exclude filtered
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa/KletKowaFil_AN_BP_DP_arenosa.vcf -R alygenomes.fasta -nt 10 --excludeFiltered -o Filtered_arenosa/KletKowaFil_AN_BP_DP_sel_arenosa.vcf 2> Filtered_arenosa/KletKowaFil_AN_BP_DP_sel_arenosa.err

#Filter 4: Allele frequency should not be 1 in both populations (fixed in both)
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_arenosa/KletKowaFil_AN_BP_DP_sel_arenosa.vcf -R alygenomes.fasta -filterName "AFfil" -filter "AF == 1.0" -o Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_arenosa.vcf 2> Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_arenosa.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_arenosa.vcf -R alygenomes.fasta -nt 10 --excludeFiltered -o Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_sel_arenosa.vcf 2> Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_sel_arenosa.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_sel_arenosa.vcf -R alygenomes.fasta -F DP -raw --allowMissingData -o Filtered_arenosa/KletKowaDP_arenosa.table 2> Filtered_arenosa/KletKowaDPtable_arenosa.err

#Filter 5: Total depth should not exceed 2.5% quantiles
grep -v "DP" Filtered_arenosa/KletKowaDP_arenosa.table > Filtered_arenosa/KletKowaDP_woheader_arenosa.table
length=$(wc -l Filtered_arenosa/KletKowaDP_woheader_arenosa.table | cut -d' ' -f1)
max=$(($length * 975 / 1000))
min=$(($length * 25 / 1000))
DPmax=$(sort -n Filtered_arenosa/KletKowaDP_woheader_arenosa.table | head -n $max | tail -n 1)
DPmin=$(sort -n Filtered_arenosa/KletKowaDP_woheader_arenosa.table | head -n $min | tail -n 1)

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_sel_arenosa.vcf -R alygenomes.fasta -filterName "DP" -filter "DP<"$DPmin".0 || DP>"$DPmax".0" -o Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_DP_arenosa.vcf 2> Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_DP_arenosa.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_DP_arenosa.vcf -R alygenomes.fasta -nt 10 --excludeFiltered -o Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_DP_sel_arenosa.vcf 2> Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_DP_sel_arenosa.err

#Split into cohorts
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_DP_sel_arenosa.vcf -R alygenomes.fasta -nt 10 -sf Kowa_arenosa.list -o Filtered_arenosa/Kowa_woGQfil_arenosa.vcf 2> Filtered_arenosa/Kowa_woGQfil_arenosa.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_DP_sel_arenosa.vcf -R alygenomes.fasta -nt 10 -sf Klet_arenosa.list -o Filtered_arenosa/Klet_woGQfil_arenosa.vcf 2> Filtered_arenosa/Klet_woGQfil_arenosa.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa/Kowa_woGQfil_arenosa.vcf -R alygenomes.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_arenosa/Kowa_woGQfil_arenosa.table 2> Filtered_arenosa/Kowa_woGQfil2_arenosa.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa/Klet_woGQfil_arenosa.vcf -R alygenomes.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_arenosa/Klet_woGQfil_arenosa.table 2> Filtered_arenosa/Klet_woGQfil2_arenosa.err
