#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan/GS_KletKowa_lyrata

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Bed_KletKowa
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

java -XX:ParallelGCThreads=2 -Xmx50g -jar ../software/GenomeAnalysisTK37.jar -T VariantsToBinaryPed -R ../alygenomes.fasta -V ../Filtered_arenosa/KletKowaFil_AN_BP_DP_AF_DP_sel_arenosa.vcf -m KletKowa.fam -bed KletKowa.bed -bim KletKowa.bim -fam KletKowa.out.fam -mgq 0 2> KletKowavcftobed.err

