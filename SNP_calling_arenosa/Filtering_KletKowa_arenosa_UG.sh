#!/bin/bash

#GATK3.7

ls UG_arenosa/*.vcf.gz | grep "Klet" | grep "vcf.gz" | sed 's/.vcf.gz//g' | sed 's/UG_arenosa\///g' > Klet_UG.list
ls UG_arenosa/*.vcf.gz | grep "Kowa" | grep "vcf.gz" | sed 's/.vcf.gz//g' | sed 's/UG_arenosa\///g' > Kowa_UG.list

#Filter 1: select bialleleic SNPs without missing genotypes
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V UG_arenosa/KletKowa_UG.vcf.gz -R alygenomes.fasta -selectType SNP -select "AN == 68" -restrictAllelesTo BIALLELIC -env -o Filtered_arenosa_UG/KletKowaFil_AN.vcf 2> Filtered_arenosa_UG/KletKowaFil_AN.err

#Filter 2: GATK best practices
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_arenosa_UG/KletKowaFil_AN.vcf -R alygenomes.fasta -filterName "BP" -filter "QD<2.0 || FS>40.0 || MQ<50.0 || MQRankSum< -2.5 || ReadPosRankSum < -4.0 || SOR>4.0" -o Filtered_arenosa_UG/KletKowaFil_AN_BP.vcf 2> Filtered_arenosa_UG/KletKowaFil_AN_BP.err

#Filter 3: Minimum depth of 2 in every sample
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_arenosa_UG/KletKowaFil_AN_BP.vcf -R alygenomes.fasta -G_filter "DP<2.0" -G_filterName "lowDP" -o Filtered_arenosa_UG/KletKowaFil_AN_DP.list 2> Filtered_arenosa_UG/KletKowaFil_AN_DP.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa_UG/KletKowaFil_AN_BP.vcf -R alygenomes.fasta -XL Filtered_arenosa_UG/KletKowaFil_AN_DP.list -o Filtered_arenosa_UG/KletKowaFil_AN_DP.vcf 2> Filtered_arenosa_UG/KletKowaFil_AN_DP2.err

grep -v "#" Filtered_arenosa_UG/KletKowaFil_AN_DP.list | grep "lowDP" >> Filtered_arenosa_UG/KletKowaFil_AN_DP.vcf 2> Filtered_arenosa_UG/KletKowaFil_lowDP.err

#Mask filter 3 to filter 2
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_arenosa_UG/KletKowaFil_AN_BP.vcf -R alygenomes.fasta --mask Filtered_arenosa_UG/KletKowaFil_AN_DP.vcf --maskName "DP" -o Filtered_arenosa_UG/KletKowaFil_AN_BP_DP.vcf 2> Filtered_arenosa_UG/KletKowaFil_AN_BP_DP.err

#Exclude filtered
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa_UG/KletKowaFil_AN_BP_DP.vcf -R alygenomes.fasta -nt 10 --excludeFiltered -o Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_sel.vcf 2> Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_sel.err

#Filter 4: Allele frequency should not be 1 in both populations (fixed in both)
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_sel.vcf -R alygenomes.fasta -filterName "AFfil" -filter "AF == 1.0" -o Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF.vcf 2> Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF.vcf -R alygenomes.fasta -nt 10 --excludeFiltered -o Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_sel.vcf 2> Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_sel.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_sel.vcf -R alygenomes.fasta -F DP -raw --allowMissingData -o Filtered_arenosa_UG/KletKowaDP.table 2> Filtered_arenosa_UG/KletKowaDPtable.err

#Filter 5: Total depth should not exceed 2.5% quantiles
grep -v "DP" Filtered_arenosa_UG/KletKowaDP.table > Filtered_arenosa_UG/KletKowaDP_woheader.table
length=$(wc -l Filtered_arenosa_UG/KletKowaDP_woheader.table | cut -d' ' -f1)
max=$(($length * 975 / 1000))
min=$(($length * 25 / 1000))
DPmax=$(sort -n Filtered_arenosa_UG/KletKowaDP_woheader.table | head -n $max | tail -n 1)
DPmin=$(sort -n Filtered_arenosa_UG/KletKowaDP_woheader.table | head -n $min | tail -n 1)

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_sel.vcf -R alygenomes.fasta -filterName "DP" -filter "DP<"$DPmin".0 || DP>"$DPmax".0" -o Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_DP.vcf 2> Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_DP.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_DP.vcf -R alygenomes.fasta -nt 10 --excludeFiltered -o Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_DP_sel.vcf 2> Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_DP_sel.err

#Split into cohorts
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_DP_sel.vcf -R alygenomes.fasta -nt 10 -sf Kowa_arenosa.list -o Filtered_arenosa_UG/Kowa_woGQfil.vcf 2> Filtered_arenosa_UG/Kowa_woGQfil.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_arenosa_UG/KletKowaFil_AN_BP_DP_AF_DP_sel.vcf -R alygenomes.fasta -nt 10 -sf Klet_arenosa.list -o Filtered_arenosa_UG/Klet_woGQfil.vcf 2> Filtered_arenosa_UG/Klet_woGQfil.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa_UG/Kowa_woGQfil.vcf -R alygenomes.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_arenosa_UG/Kowa_UG.table 2> Filtered_arenosa_UG/Kowa_woGQfil2.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_arenosa_UG/Klet_woGQfil.vcf -R alygenomes.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_arenosa_UG/Klet_UG.table 2> Filtered_arenosa_UG/Klet_woGQfil2.err
