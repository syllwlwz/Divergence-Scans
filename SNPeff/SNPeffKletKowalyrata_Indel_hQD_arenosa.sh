#!/bin/bash

java -Xmx50g -jar software/snpEff/snpEff.jar LyV2 Indel_filtered_HC_new_arenosa/KletKowa_INDEL_BP_DP_AF_DP_sel_hQD.vcf > SNPeff/KletKowa_INDEL_hQD.arenosa.ann.vcf 2> SNPeff/KletKowa_INDEL_hQD.arenosa.ann.vcf.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V SNPeff/KletKowa_INDEL_hQD.arenosa.ann.vcf -R alygenomes.fasta -F CHROM -F POS -F ANN -raw --allowMissingData -o SNPeff/KletKowa_INDEL_hQD_ann.arenosa.table 2> SNPeff/KletKowa_INDEL_hQD_ann.arenosa.table.err
