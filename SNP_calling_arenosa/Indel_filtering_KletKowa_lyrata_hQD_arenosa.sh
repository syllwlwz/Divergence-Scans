#!/bin/bash

#GATK3.7

#ls HC1_lyrata/*.vcf.gz | grep "Klet" | grep "vcf.gz" | sed 's/.vcf.gz//g' | sed 's/HC1_lyrata\///g' > Klet.list
#ls HC1_lyrata/*.vcf.gz | grep "Kowa" | grep "vcf.gz" | sed 's/.vcf.gz//g' | sed 's/HC1_lyrata\///g' > Kowa.list

#Filter 1: select bialleleic SNPs without missing genotypes
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V HC2_arenosa/KletKowa_Indel.vcf -R alygenomes.fasta -selectType INDEL -restrictAllelesTo BIALLELIC -select "AN == 68" -o Indel_filtered_HC_new_arenosa/KletKowa_INDEL_hQD.vcf 2> Indel_filtered_HC_new_arenosa/KletKowa_INDEL_hQD.err

#Filter 2: GATK best practices
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Indel_filtered_HC_new_arenosa/KletKowa_INDEL_hQD.vcf -R alygenomes.fasta -filterName "BP" -filter "QD<20.0 || FS>40.0 || ReadPosRankSum < -4.0 || SOR>4.0" -o Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_hQD.vcf 2> Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_hQD.err

#Filter 3: Minimum depth of 2 in every sample
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_hQD.vcf -R alygenomes.fasta -G_filter "DP<2.0" -G_filterName "lowDP" -o Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_hQD.list 2> Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_hQD.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_hQD.vcf -R alygenomes.fasta -XL Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_hQD.list -o Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_hQD.vcf 2> Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP2_hQD.err

grep -v "#" Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_hQD.list | grep "lowDP" >> Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_hQD.vcf 2> Indel_filtered_HC_new_arenosa/KletKowa_lowDP_hQD.err

#Mask filter 3 to filter 2
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_hQD.vcf -R alygenomes.fasta --mask Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_hQD.vcf --maskName "DP" -o Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_hQD.vcf 2> Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_hQD.err

#Exclude filtered
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_hQD.vcf -R alygenomes.fasta -nt 10 --excludeFiltered -o Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_sel_hQD.vcf 2> Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_sel_hQD.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_sel_hQD.vcf -R alygenomes.fasta -F DP -raw --allowMissingData -o Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_hQD.table 2> Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DPtable_hQD.err

#Filter 5: Total depth should not exceed 2.5% quantiles
grep -v "DP" Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_hQD.table > Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_woheader_hQD.table
length=$(wc -l Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_woheader_hQD.table | cut -d' ' -f1)
max=$(($length * 975 / 1000))
min=$(($length * 25 / 1000))
DPmax=$(sort -n Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_woheader_hQD.table | head -n $max | tail -n 1)
DPmin=$(sort -n Indel_filtered_HC_new_arenosa/KletKowa_INDEL_DP_woheader_hQD.table | head -n $min | tail -n 1)

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_sel_hQD.vcf -R alygenomes.fasta -filterName "DP" -filter "DP<"$DPmin".0 || DP>"$DPmax".0" -o Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_AF_DP_hQD.vcf 2> Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_AF
_DP_hQD.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_AF_DP_hQD.vcf -R alygenomes.fasta -nt 10 --excludeFiltered -o Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_AF_DP_sel_hQD.vcf 2> Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_AF_DP_sel_hQD.err

#Split into cohorts
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_AF_DP_sel_hQD.vcf -R alygenomes.fasta -nt 10 -sf Kowa_arenosa.list -o Indel_filtered_HC_new_arenosa/Kowa_INDEL_hQD.vcf 2> Indel_filtered_HC_new_arenosa/Kowa_INDEL_hQD.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_AF_DP_sel_hQD.vcf -R alygenomes.fasta -nt 10 -sf Klet_arenosa.list -o Indel_filtered_HC_new_arenosa/Klet_INDEL_hQD.vcf 2> Indel_filtered_HC_new_arenosa/Klet_INDEL_hQD.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Indel_filtered_HC_new_arenosa/Kowa_INDEL_hQD.vcf -R alygenomes.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Indel_filtered_HC_new_arenosa/Kowa_INDEL_hQD.table 2> Indel_filtered_HC_new_arenosa/Kowa_INDEL2_hQD.err

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Indel_filtered_HC_new_arenosa/Klet_INDEL_hQD.vcf -R alygenomes.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Indel_filtered_HC_new_arenosa/Klet_INDEL_hQD.table 2> Indel_filtered_HC_new_arenosa/Klet_INDEL2_hQD.err

