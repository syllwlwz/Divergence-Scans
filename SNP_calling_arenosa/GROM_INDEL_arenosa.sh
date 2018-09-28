#!/bin/bash

#Select INDELs


grep "SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP\|SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP"  GROM_analysis/Klet_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Klet_indels_arenosa.vcf

grep "SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP\|SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP"  GROM_analysis/Kowa_arenosa.predictions.vcf | grep -v "#" > GROM_analysis/Kowa_indels_arenosa.vcf

