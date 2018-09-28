#!/bin/bash

java -Xmx50g -jar software/snpEff/snpEff.jar LyV2 Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_DP_sel_arenosa.vcf -t > SNPeff_scans/KletKowa_arenosa.ann.vcf 2> SNPeff_scans/KletKowa_arenosa.ann.vcf.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V SNPeff_scans/KletKowa_arenosa.ann.vcf -R alygenomes.fasta -F CHROM -F POS -F ANN -raw --allowMissingData -o SNPeff_scans/KletKowa_arenosa.ann.table 2> SNPeff_scans/KletKowa_arenosa.ann.table.err
