# Divergence-Scans
Genome scans and high effect indel and SNP selection

example for Klet/ Klet and Kowa in arenosa
replace arenosa samples with halleri samples and Klet for other population names as well as parameters as necessary
file paths will have to be changed to actual locations of files

SNP and indel calling
1. Cutadapt_arenosa.sh
2. Alignment_lyrata_arenosa.sh
3. Deduplication_arenosa.sh
4. Namefix_arenosa.sh
5. Realignment_arenosa.sh
HaplotypeCaller: 
6.a) HC1_arenosa.sh b) HC2_KletKowa_arenosa_list.sh
7. For SNPs:Filtering_arenosa.sh
8. For Indels: Indel_filtering_KletKowa_arenosa_test.sh
8. Filtering_plots.R
9. for SNPS: Filtering_KletKowa_arenosa.sh
9. For Indels: Indel_filtering_KletKowa_lyrata_hQD_arenosa.sh

UnifiedGenotyper:
6. UG_arenosa.sh
9. Filtering_KletKowa_arenosa_UG.sh

GROM:
6. GROM_arenosa.sh
7. GROM_INDEL_arenosa.sh

Indel analysis
1. SNPeffKletKowalyrata_Indel_GROM_arenosa.sh
2. SNPeffKletKowalyrata_Indel_hQD_arenosa.sh
3. Indel_analysis_KK_hQD_arenosa.R
4. Coverage: a) Hist_coverage_arenosa.sh b) Hist_coverage_arenosa_exons.sh
5. GROM_analysis_arenosa.R (for all populations)

Genome scans
1. Post_HC_analysis_KletKowa_woGQfil_arenosa.R
2.  DD_KletKowalyrata_woGQfil.R
	Dxy_KletKowa_arenosa.R
	FayandWusH_KletKowa_lyrata_mean_arenosa.R
	FLK_KletKowalyrata_woGQfil_min_arenosa_alt.R
	FST_KletKowalyrata_arenosa_corrected.R
	NielsenKletKowalyrata_arenosa.R
	Tajimas_D_KletKowa_lyrata_mean_arenosa.R
	VarLD: a) KletKowa_GT_arenosa.sh b) VarLD_KletKowa_arenosa.r c) KletKowa_VarLD.sh d) KletKowa_VarLD_25SNPwindows_arenosa.R
	Flk: a) KletKowa_bed_arenosa.sh b) FLK_KletKowalyrata_woGQfil_min_arenosa_alt.R
	SweeD: SweeD_Klet.sh (Afterwards to split files for each scaffold :csplit -f Klet_Sweed SweeD_Report.Klet_SweeD "///+1" "{*}" and head -n -2 Klet_Sweed01 > Klet_Sweed_scaffold1.table and combine with paste in unix)
	AF_KletKowa_lyrata.R

Genome scan analysis
1. SNPeffKletKowalyrata_arenosa.sh
2. Genome_scan_analysis_with_promoters_KletKowa_1percentoutliers_lyrata_TDmean_Flkmin_arenosa.R
	for allele frequency plots: AF_Plots_Klet-Kowa_arenosa.R
	for metrics plots: Metrics_plots_KletKowa_arenosa.R
3. Additionalinfo_GSfile_01percent_KK_arenosa.R
4. Comparison_arenosa.R to find overlap between population pairs
5. Comparisonbetweenspecies.R to find overlap between all population pairs, to add Araport info and to reformat tables (additional formatting done in Excel e.g. renaming of columns)
