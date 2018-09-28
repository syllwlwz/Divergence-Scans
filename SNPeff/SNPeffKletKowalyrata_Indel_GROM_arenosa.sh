#!/bin/bash

java -Xmx50g -jar software/snpEff/snpEff.jar LyV2 GROM_analysis/Klet_indels_arenosa.vcf > SNPeff/KletKowa_INDEL_GROM.arenosa.ann.vcf 2> SNPeff/KletKowa_INDEL_GROM.arenosa.ann.vcf.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V SNPeff/KletKowa_INDEL_GROM.arenosa.ann.vcf -R alygenomes.fasta -F CHROM -F POS -F ANN -raw --allowMissingData -o SNPeff/KletKowa_INDEL_GROM.arenosa.ann.table 2> SNPeff/KletKowa_INDEL_GROM.arenosa.ann.table.err


java -Xmx50g -jar software/snpEff/snpEff.jar LyV2 GROM_analysis/Kowa_indels_arenosa.vcf > SNPeff/KletKowa_INDEL_GROM.Kowa.arenosa.ann.vcf 2> SNPeff/KletKowa_INDEL_GROM.Kowa.arenosa.ann.vcf.err

#Extract tables
java -XX:ParallelGCThreads=2 -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V SNPeff/KletKowa_INDEL_GROM.Kowa.arenosa.ann.vcf -R alygenomes.fasta -F CHROM -F POS -F ANN -raw --allowMissingData -o SNPeff/KletKowa_INDEL_GROM.Kowa.arenosa.ann.table 2> SNPeff/KletKowa_INDEL_GROM_Kowa.table.err


