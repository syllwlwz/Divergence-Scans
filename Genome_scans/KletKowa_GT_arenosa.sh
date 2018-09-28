#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan/GS_KletKowa_lyrata

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N GT_KletKowa
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

java -XX:ParallelGCThreads=2 -Xmx50g -jar ../software/GenomeAnalysisTK37.jar -T VariantsToTable -V ../Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_DP_sel_arenosa.vcf -R ../alygenomes.fasta -F CHROM -F POS -GF GT -raw --allowMissingData -o KletKowa_GT.table 2> KletKowa_GT.err
