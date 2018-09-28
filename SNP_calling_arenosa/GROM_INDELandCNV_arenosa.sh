#!/bin/bash

#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N GROM_Filtering_arenosa
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#Select INDELs and CNVs

#grep "SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP\|SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP"  GROM_analysis/Mias_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Mias_indels_arenosa.vcf
#grep "SD:Z:CN:CS"  GROM_analysis/Mias_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Mias_cnvs_arenosa.vcf

grep "SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP\|SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP"  GROM_analysis/Zapa_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Zapa_indels_arenosa.vcf
#grep "SD:Z:CN:CS"  GROM_analysis/Zapa_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Zapa_cnvs_arenosa.vcf

grep "SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP\|SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP"  GROM_analysis/Klet_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Klet_indels_arenosa.vcf
#grep "SD:Z:CN:CS"  GROM_analysis/Klet_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Klet_cnvs_arenosa.vcf

grep "SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP\|SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP"  GROM_analysis/Kowa_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Kowa_indels_arenosa.vcf
#grep "SD:Z:CN:CS"  GROM_analysis/Kowa_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Kowa_cnvs_arenosa.vcf

grep "SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP\|SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP"  GROM_analysis/Wulm_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Wulm_indels_arenosa.vcf
#grep "SD:Z:CN:CS"  GROM_analysis/Wulm_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Wulm_cnvs_arenosa.vcf

grep "SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP\|SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP"  GROM_analysis/Chok_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Chok_indels_arenosa.vcf
#grep "SD:Z:CN:CS"  GROM_analysis/Chok_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Chok_cnvs_arenosa.vcf


