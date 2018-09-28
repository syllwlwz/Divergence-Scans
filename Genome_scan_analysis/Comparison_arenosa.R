require(xlsx)
options(java.parameters = "-Xmx16g")

#MiasZapa directory
MZ_FST<-read.csv("FST_MiasZapa_01percent_arenosa.csv",header=T)
MZ_DXY<-read.csv("DXY_MiasZapa_01percent_arenosa.csv",header=T)
MZ_AFDabs<-read.csv("AFDabs_MiasZapa_01percent_arenosa.csv",header=T)
MZ_NIELSEN<-read.csv("NIELSEN_MiasZapa_01percent_arenosa.csv",header=T)
MZ_FLK<-read.csv("FLK_MiasZapa_01percent_arenosa.csv",header=T)
MZ_DD<-read.table("DD_MiasZapa_01percent_arenosa.csv",header=T,sep=";")
MZ_HighSNP<-read.table("GenesHIGHeffect_above40percentAFdiff_withUG_MZ_arenosa.txt", sep="\t",header=T)
MZ_HighSNP<-MZ_HighSNP[MZ_HighSNP$UG_AF_diff!="NO",]

#change to GROM_analysis directory
MZ_INDEL<-read.table("Indel_MZ_HC_GROM_arenosa.table",sep="\t",header=T)

MZ_FSTr<-cbind(MZ_FST,rep("FST",nrow(MZ_FST)))
MZ_DXYr<-cbind(MZ_DXY,rep("DXY",nrow(MZ_DXY)))
MZ_AFDabsr<-cbind(MZ_AFDabs,rep("AFDabs",nrow(MZ_AFDabs)))
MZ_NIELSENr<-cbind(MZ_NIELSEN,rep("NIELSEN",nrow(MZ_NIELSEN)))
MZ_FLKr<-cbind(MZ_FLK,rep("FLK",nrow(MZ_FLK)))
MZ_DDr<-cbind(MZ_DD,rep("DD",nrow(MZ_DD)))

names(MZ_FSTr)[69]<-"Outlier_metric"
names(MZ_DXYr)[69]<-"Outlier_metric"
names(MZ_AFDabsr)[69]<-"Outlier_metric"
names(MZ_NIELSENr)[69]<-"Outlier_metric"
names(MZ_FLKr)[69]<-"Outlier_metric"
names(MZ_DDr)[69]<-"Outlier_metric"

names(MZ_FSTr)[50]<-"Rank_in_outlier_metric"
names(MZ_DXYr)[50]<-"Rank_in_outlier_metric"
names(MZ_AFDabsr)[50]<-"Rank_in_outlier_metric"
names(MZ_NIELSENr)[50]<-"Rank_in_outlier_metric"
names(MZ_FLKr)[50]<-"Rank_in_outlier_metric"
names(MZ_DDr)[50]<-"Rank_in_outlier_metric"

names(MZ_FLKr)[18]<-"Flk_pvalue"
MZ<-rbind(MZ_FSTr,MZ_DXYr,MZ_AFDabsr,MZ_NIELSENr,MZ_FLKr,MZ_DDr)

MZ[,9:11]<-apply(MZ[,9:11],2,function(x) as.numeric(as.character(x)))
MZ[,13:25]<-apply(MZ[,13:25],2,function(x) as.numeric(as.character(x)))
MZ[,38:55]<-apply(MZ[,38:55],2,function(x) as.numeric(as.character(x)))

MZFiltered<-MZ[MZ$AF_Plot!="D"&MZ$Metrics_Plot!="C"&(MZ$TD_diff_nearest_Mias_Zapa_stand<0|MZ$FWH_diff_nearest_Mias_Zapa_stand<0|MZ$Sweed_Mias_nearest_Likelihood>MZ$Sweed_Zapa_nearest_Likelihood)&(MZ$TajimasD_Mias_nearest_stand<0|MZ$FayandWusH_Mias_nearest_stand<0|MZ$Sweed_Mias_nearest_Likelihood>0),]
MZFiltered<-MZFiltered[!duplicated(MZFiltered[,c(1,6:10,13:25,38:55)]),]
#write.table(MZFiltered,"TestMZ.table",sep="\t",row.names=F)

FST<-read.csv("Summary_genes_1percent_MiasZapa_FST_arenosa.csv",header=T,sep=";")
DXY<-read.csv("Summary_genes_1percent_MiasZapa_DXY_arenosa.csv",header=T,sep=";")
AFDabs<-read.csv("Summary_genes_1percent_MiasZapa_AFDabs_arenosa.csv",header=T,sep=";")
NIELSEN<-read.csv("Summary_genes_1percent_MiasZapa_NIELSEN_arenosa.csv",header=T,sep=";")
FLK<-read.csv("Summary_genes_1percent_MiasZapa_FLK_arenosa.csv",header=T,sep=";")
DD<-read.csv("Summary_genes_1percent_MiasZapa_DD_arenosa.csv",header=T,sep=";")

FST<-FST[,1:30]
DXY<-DXY[,1:30]
AFDabs<-AFDabs[,1:30]
NIELSEN<-NIELSEN[,1:30]
FLK<-FLK[,1:30]
DD<-DD[,1:30]

FST<-cbind(FST,rep("FST",nrow(FST)))
DXY<-cbind(DXY,rep("DXY",nrow(DXY)))
AFDabs<-cbind(AFDabs,rep("AFDabs",nrow(AFDabs)))
NIELSEN<-cbind(NIELSEN,rep("NIELSEN",nrow(NIELSEN)))
FLK<-cbind(FLK,rep("FLK",nrow(FLK)))
DD<-cbind(DD,rep("DD",nrow(DD)))

names(FST)[31]<-"Metric"
names(DXY)[31]<-"Metric"
names(AFDabs)[31]<-"Metric"
names(NIELSEN)[31]<-"Metric"
names(FLK)[31]<-"Metric"
names(DD)[31]<-"Metric"

FST$Rank<-rank(-as.numeric(as.character(FST$Fst)),ties.method="min")
NIELSEN$Rank<-rank(-as.numeric(as.character(NIELSEN$Nielsen)),ties.method="min")
DXY$Rank<-rank(-as.numeric(as.character(DXY$Dxy)),ties.method="min")
DD$Rank<-rank(as.numeric(as.character(DD$DD)),ties.method="min")
AFDabs$Rank<-rank(-as.numeric(as.character(AFDabs$AFDabs)),ties.method="min")
FLK$Rank<-rank(-as.numeric(as.character(FLK$Flk)),ties.method="min")

all<-rbind(FST,DXY,AFDabs,NIELSEN,FLK,DD)
all<-all[!duplicated(all[,c(1,7,8,9,10,32)]),]
#MZFiltered<-MZFiltered[!duplicated(MZFiltered[,c(1,7,8,9,10,74)]),]

test<-ftable(Gene~Metric,data=all)
Numberofmetrics1percent<-colSums(test!=0)
Genesnummet<-c(unlist(attr(test,"col.vars"),use.names=F))
Numberofmetrics1percentgenes<-as.data.frame(cbind(Genesnummet,Numberofmetrics1percent))
names(Numberofmetrics1percentgenes)<-c("Gene","Numberofmetrics1percent")
Numberofmetrics1percentgenes[,1]<-as.character(Numberofmetrics1percentgenes[,1])
Numberofmetrics1percentgenes[,2]<-as.numeric(as.character(Numberofmetrics1percentgenes[,2]))
Numberofmetrics1percentgenes[Numberofmetrics1percentgenes$Gene=="",]<-NA

MZFiltered1percentmerge<-merge(MZFiltered,Numberofmetrics1percentgenes,by="Gene",all.x=T)
MZFiltered1percentmerge<-MZFiltered1percentmerge[!duplicated(MZFiltered1percentmerge[,c(1,6:10,13:30,34,38:55,69,70)]),]

write.table(MZFiltered1percentmerge,"MiasZapa_Filtered_w1percentoverlap_incl.table",sep="\t",row.names=F)

length(unique(MZFiltered1percentmerge$Gene))
#163



########################################
#KletKowa#
########################################
#KletKowa directory
KK_FST<-read.csv("FST_KletKowa_01percent_arenosa.csv",header=T)
KK_DXY<-read.csv("DXY_KletKowa_01percent_arenosa.csv",header=T)
KK_AFDabs<-read.csv("AFDabs_KletKowa_01percent_arenosa.csv",header=T)
KK_NIELSEN<-read.csv("NIELSEN_KletKowa_01percent_arenosa.csv",header=T)
KK_FLK<-read.csv("FLK_KletKowa_01percent_arenosa.csv",header=T)
KK_DD<-read.csv("DD_KletKowa_01percent_arenosa.csv",header=T)
KK_HighSNP<-read.table("GenesHIGHeffect_above40percentAFdiff_withUG_KK_arenosa.txt", sep="\t",header=T)
KK_HighSNP<-KK_HighSNP[KK_HighSNP$UG_AF_diff!="NO",]

#change to GROM_analysis directory
KK_INDEL<-read.table("Indel_KK_HC_GROM_arenosa.table",sep="\t",header=T)

KK_FSTr<-cbind(KK_FST,rep("FST",nrow(KK_FST)))
KK_DXYr<-cbind(KK_DXY,rep("DXY",nrow(KK_DXY)))
KK_AFDabsr<-cbind(KK_AFDabs,rep("AFDabs",nrow(KK_AFDabs)))
KK_NIELSENr<-cbind(KK_NIELSEN,rep("NIELSEN",nrow(KK_NIELSEN)))
KK_FLKr<-cbind(KK_FLK,rep("FLK",nrow(KK_FLK)))
KK_DDr<-cbind(KK_DD,rep("DD",nrow(KK_DD)))

names(KK_FSTr)[69]<-"Outlier_metric"
names(KK_DXYr)[69]<-"Outlier_metric"
names(KK_AFDabsr)[69]<-"Outlier_metric"
names(KK_NIELSENr)[69]<-"Outlier_metric"
names(KK_FLKr)[69]<-"Outlier_metric"
names(KK_DDr)[69]<-"Outlier_metric"

names(KK_FSTr)[50]<-"Rank_in_outlier_metric"
names(KK_DXYr)[50]<-"Rank_in_outlier_metric"
names(KK_AFDabsr)[50]<-"Rank_in_outlier_metric"
names(KK_NIELSENr)[50]<-"Rank_in_outlier_metric"
names(KK_FLKr)[50]<-"Rank_in_outlier_metric"
names(KK_DDr)[50]<-"Rank_in_outlier_metric"

names(KK_FLKr)[18]<-"Flk_pvalue"
KK<-rbind(KK_FSTr,KK_DXYr,KK_AFDabsr,KK_NIELSENr,KK_FLKr,KK_DDr)

KK[,9:11]<-apply(KK[,9:11],2,function(x) as.numeric(as.character(x)))
KK[,13:25]<-apply(KK[,13:25],2,function(x) as.numeric(as.character(x)))
KK[,38:55]<-apply(KK[,38:55],2,function(x) as.numeric(as.character(x)))

KKFiltered<-KK[KK$AF_Plot!="D"&KK$Metrics_Plot!="C"&(KK$TD_diff_nearest_Klet_Kowa_stand<0|KK$FWH_diff_nearest_Klet_Kowa_stand<0|KK$Sweed_Klet_nearest_Likelihood>KK$Sweed_Kowa_nearest_Likelihood)&(KK$TajimasD_Klet_nearest_stand<0|KK$FayandWusH_Klet_nearest_stand<0|KK$Sweed_Klet_nearest_Likelihood>0),]
KKFiltered<-KKFiltered[!duplicated(KKFiltered[,c(1,6:10,13:25,38:55)]),]
#write.table(KKFiltered,"TestKK.table",sep="\t",row.names=F)

FST<-read.csv("Summary_genes_1percent_KletKowa_FST_arenosa.csv",header=T,sep=";")
DXY<-read.csv("Summary_genes_1percent_KletKowa_DXY_arenosa.csv",header=T,sep=";")
AFDabs<-read.csv("Summary_genes_1percent_KletKowa_AFDabs_arenosa.csv",header=T,sep=";")
NIELSEN<-read.csv("Summary_genes_1percent_KletKowa_NIELSEN_arenosa.csv",header=T,sep=";")
FLK<-read.csv("Summary_genes_1percent_KletKowa_FLK_arenosa.csv",header=T,sep=";")
DD<-read.csv("Summary_genes_1percent_KletKowa_DD_arenosa.csv",header=T,sep=";")

FST<-FST[,1:30]
DXY<-DXY[,1:30]
AFDabs<-AFDabs[,1:30]
NIELSEN<-NIELSEN[,1:30]
FLK<-FLK[,1:30]
DD<-DD[,1:30]

FST<-cbind(FST,rep("FST",nrow(FST)))
DXY<-cbind(DXY,rep("DXY",nrow(DXY)))
AFDabs<-cbind(AFDabs,rep("AFDabs",nrow(AFDabs)))
NIELSEN<-cbind(NIELSEN,rep("NIELSEN",nrow(NIELSEN)))
FLK<-cbind(FLK,rep("FLK",nrow(FLK)))
DD<-cbind(DD,rep("DD",nrow(DD)))

names(FST)[31]<-"Metric"
names(DXY)[31]<-"Metric"
names(AFDabs)[31]<-"Metric"
names(NIELSEN)[31]<-"Metric"
names(FLK)[31]<-"Metric"
names(DD)[31]<-"Metric"

FST$Rank<-rank(-as.numeric(as.character(FST$Fst)),ties.method="min")
NIELSEN$Rank<-rank(-as.numeric(as.character(NIELSEN$Nielsen)),ties.method="min")
DXY$Rank<-rank(-as.numeric(as.character(DXY$Dxy)),ties.method="min")
DD$Rank<-rank(as.numeric(as.character(DD$DD)),ties.method="min")
AFDabs$Rank<-rank(-as.numeric(as.character(AFDabs$AFDabs)),ties.method="min")
FLK$Rank<-rank(-as.numeric(as.character(FLK$Flk)),ties.method="min")

all<-rbind(FST,DXY,AFDabs,NIELSEN,FLK,DD)
all<-all[!duplicated(all[,c(1,7,8,9,10,32)]),]
#KKFiltered<-KKFiltered[!duplicated(KKFiltered[,c(1,7,8,9,10,74)]),]

test<-ftable(Gene~Metric,data=all)
Numberofmetrics1percent<-colSums(test!=0)
Genesnummet<-c(unlist(attr(test,"col.vars"),use.names=F))
Numberofmetrics1percentgenes<-as.data.frame(cbind(Genesnummet,Numberofmetrics1percent))
names(Numberofmetrics1percentgenes)<-c("Gene","Numberofmetrics1percent")
Numberofmetrics1percentgenes[,1]<-as.character(Numberofmetrics1percentgenes[,1])
Numberofmetrics1percentgenes[,2]<-as.numeric(as.character(Numberofmetrics1percentgenes[,2]))
Numberofmetrics1percentgenes[Numberofmetrics1percentgenes$Gene=="",]<-NA

KKFiltered1percentmerge<-merge(KKFiltered,Numberofmetrics1percentgenes,by="Gene",all.x=T)
KKFiltered1percentmerge<-KKFiltered1percentmerge[!duplicated(KKFiltered1percentmerge[,c(1,6:10,13:30,34,38:55,69,70)]),]

write.table(KKFiltered1percentmerge,"KletKowa_Filtered_w1percentoverlap_incl.table",sep="\t",row.names=F)

length(unique(KKFiltered1percentmerge$Gene))
#173


HM_list<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",2,header=T)
KKFiltered1percentmerge_HMnew<-merge(KKFiltered1percentmerge[,-c(56:68)],HM_list,by.x="Gene",by.y="AGI_code",all.x=T)

#write.table(KKFiltered1percentmerge,"KletKowa_Filtered_w1percentoverlap_incl_newHM.table",sep="\t",row.names=F)

MZFiltered1percentmerge_HMnew<-merge(MZFiltered1percentmerge[,-c(56:68)],HM_list,by.x="Gene",by.y="AGI_code",all.x=T)

#write.table(MZFiltered1percentmerge,"MiasZapa_Filtered_w1percentoverlap_incl_newHM.table",sep="\t",row.names=F)

#add other lists

#Grom_analysis directory
family_info<-read.table("gene_families_araport_sel.txt",sep="\t",header=T,quote="",fill=T)
family_info$Genomic_Locus_Tag<-toupper(as.character(family_info$Genomic_Locus_Tag))

MZFiltered1percentmerge_HMnew_fam<-merge(MZFiltered1percentmerge_HMnew,family_info,by.x="Gene",by.y="Genomic_Locus_Tag",all.x=T)
KKFiltered1percentmerge_HMnew_fam<-merge(KKFiltered1percentmerge_HMnew,family_info,by.x="Gene",by.y="Genomic_Locus_Tag",all.x=T)


###########################################################################
#Comparison#
###########################################################################
require(GenomicRanges)

KK_HighSNP_fixed<-KK_HighSNP[KK_HighSNP$AF_Klet<=0.1|KK_HighSNP$AF_Klet>=0.9,]
MZ_HighSNP_fixed<-MZ_HighSNP[MZ_HighSNP$AF_Mias<=0.1|MZ_HighSNP$AF_Mias>=0.9,]

Klet_cn<-read.table("Klet_CNVs_overlap_strict.table",header=T,sep="\t")
Mias_cn<-read.table("Mias_CNVs_overlap_strict.table",header=T,sep="\t")

KK_classical_overlap<-KKFiltered1percentmerge_HMnew_fam[,1:7]
KK_classical_overlap$Metric<-KKFiltered1percentmerge_HMnew_fam$Outlier_metric
KK_Indel_overlap<-KK_INDEL[,1:7]
KK_Indel_overlap$Metric<-rep("KK_INDEL",nrow(KK_Indel_overlap))
KK_higheffectSNP_overlap<-KK_HighSNP_fixed[,1:7]
KK_higheffectSNP_overlap$Metric<-rep("KK_higheffectSNP",nrow(KK_higheffectSNP_overlap))

KK_overlap<-rbind(KK_classical_overlap,KK_higheffectSNP_overlap,KK_Indel_overlap)
KK_overlap<-KK_overlap[!duplicated(KK_overlap),]


MZ_classical_overlap<-MZFiltered1percentmerge_HMnew_fam[,1:7]
MZ_classical_overlap$Metric<-MZFiltered1percentmerge_HMnew_fam$Outlier_metric
MZ_Indel_overlap<-MZ_INDEL[,1:7]
MZ_Indel_overlap$Metric<-rep("MZ_INDEL",nrow(MZ_Indel_overlap))
MZ_higheffectSNP_overlap<-MZ_HighSNP_fixed[,1:7]
MZ_higheffectSNP_overlap$Metric<-rep("MZ_higheffectSNP",nrow(MZ_higheffectSNP_overlap))

MZ_overlap<-rbind(MZ_classical_overlap,MZ_higheffectSNP_overlap,MZ_Indel_overlap)
MZ_overlap<-MZ_overlap[!duplicated(MZ_overlap),]

KK_overlap_all<-KK_overlap
MZ_overlap_all<-MZ_overlap

KK_overlap_all_fam<-merge(KK_overlap_all,family_info,by.x="Gene",by.y="Genomic_Locus_Tag")
MZ_overlap_all_fam<-merge(MZ_overlap_all,family_info,by.x="Gene",by.y="Genomic_Locus_Tag")

KK_MZ_overlap_all_fam<-merge(KK_overlap_all_fam,MZ_overlap_all_fam,by="Gene_Family")

KK_MZ_overlap_all_fam<-KK_MZ_overlap_all_fam[!duplicated(KK_MZ_overlap_all_fam),]
names(KK_MZ_overlap_all_fam)[9]<-"Metric_KK"
names(KK_MZ_overlap_all_fam)[18]<-"Metric_MZ"

Methomlist<-HM_list
Methomlist$AGI_Number<-toupper(Methomlist$AGI_code)

overlap_all_fam_renamed_KK<-overlap_all_fam[,1:10]
overlap_all_fam_renamed_MZ<-overlap_all_fam[,c(1,11:19)]
overlap_all_fam_renamed_WH<-overlap_all_fam[,c(1,20:28)]

names(overlap_all_fam_renamed_KK)<-c("Gene_Family","Gene","Type","Short_description","Curator_summary","Computational_description","OrthoGroup","Lyr_Gene_Ahal_Gene","Metric","Sub_Family")
overlap_all_fam_renamed_KK$Population_pair<-rep("KK",nrow(overlap_all_fam_renamed_KK))
names(overlap_all_fam_renamed_MZ)<-c("Gene_Family","Gene","Type","Short_description","Curator_summary","Computational_description","OrthoGroup","Lyr_Gene_Ahal_Gene","Metric","Sub_Family")
overlap_all_fam_renamed_MZ$Population_pair<-rep("MZ",nrow(overlap_all_fam_renamed_MZ))
names(overlap_all_fam_renamed_WH)<-c("Gene_Family","Gene","Type","Short_description","Curator_summary","Computational_description","OrthoGroup","Lyr_Gene_Ahal_Gene","Metric","Sub_Family")
overlap_all_fam_renamed_WH$Population_pair<-rep("WH",nrow(overlap_all_fam_renamed_WH))

all<-rbind(overlap_all_fam_renamed_KK,overlap_all_fam_renamed_MZ,overlap_all_fam_renamed_WH)
all_overlap<-unique(all)
all_overlap_metlist<-merge(all_overlap,Methomlist,all.x=T,by.x="Gene",by.y="AGI_Number")
write.table(all_overlap_metlist,"Overlap_all_fam_morphed.table",sep="\t",row.names=F)

#Overlap thaliana genes

KK_overlap_all<-na.exclude(KK_overlap_all[!KK_overlap_all$Gene=="",])
MZ_overlap_all<-na.exclude(MZ_overlap_all[!MZ_overlap_all$Gene=="",])

KK_MZ_overlap<-merge(KK_overlap_all,MZ_overlap_all,by="Gene")

KK_MZ_overlap<-KK_MZ_overlap[!duplicated(KK_MZ_overlap),]
names(KK_MZ_overlap)[8]<-"Metric_KK"
names(KK_MZ_overlap)[15]<-"Metric_MZ"

KK_MZ_overlap_metlist<-merge(KK_MZ_overlap,Methomlist,all.x=T,by.x="Gene",by.y="AGI_code")

write.table(unique(KK_MZ_overlap_metlist),"KK_MZ_overlap.table",sep="\t",row.names=F)

#lyrata genes or halleri genesoverlapping

MZ_overlap_all<-MZ_overlap
KK_overlap_all<-KK_overlap

KK_MZ_overlap_lyrgene<-unique(merge(KK_overlap_all,MZ_overlap_all[,7:8],by="Lyr_Gene"))

#Methomlist<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",2,header=T)
#Methomlist$AGI_Number<-toupper(Methomlist$AGI_code)

KK_MZ_overlap_lyrgene_metlist<-merge(KK_MZ_overlap_lyrgene,Methomlist,all.x=T,by.x="Gene",by.y="AGI_Number")

write.table(KK_MZ_overlap_lyrgene_metlist,"KK_MZ_overlap_lyrgene.table",sep="\t",row.names=F)

options(java.parameters = "-Xmx1024m")

#GROM_analysis
Mias_coverage<-read.table("Mias_arenosa_genes_exons.table",sep="\t",header=T)
Zapa_coverage<-read.table("Zapa_arenosa_genes_exons.table",sep="\t",header=T)
Klet_coverage<-read.table("Klet_arenosa_genes_exons.table",sep="\t",header=T)
Kowa_coverage<-read.table("Kowa_arenosa_genes_exons.table",sep="\t",header=T)

#Comparison
MZFiltered1percentmerge_cov1<-merge(MZFiltered1percentmerge,Mias_coverage[,c(1,10:13)],by.x="Lyr_Gene",by.y="ID")
MZFiltered1percentmerge_cov2<-merge(MZFiltered1percentmerge_cov1,Zapa_coverage[,c(1,10:13)],by.x="Lyr_Gene",by.y="ID")
names(MZFiltered1percentmerge_cov2)[71:78]<-c("Mean_coverage_gene_Mias","Median_coverage_gene_Mias","Mean_coverage_exons_Mias","Median_coverage_exons_Mias","Mean_coverage_gene_Zapa","Median_coverage_gene_Zapa","Mean_coverage_exons_Zapa","Median_coverage_exons_Zapa")
KKFiltered1percentmerge_cov1<-merge(KKFiltered1percentmerge,Klet_coverage[,c(1,10:13)],by.x="Lyr_Gene",by.y="ID")
KKFiltered1percentmerge_cov2<-merge(KKFiltered1percentmerge_cov1,Kowa_coverage[,c(1,10:13)],by.x="Lyr_Gene",by.y="ID")
names(KKFiltered1percentmerge_cov2)[71:78]<-c("Mean_coverage_gene_Klet","Median_coverage_gene_Klet","Mean_coverage_exons_Klet","Median_coverage_exons_Klet","Mean_coverage_gene_Kowa","Median_coverage_gene_Kowa","Mean_coverage_exons_Kowa","Median_coverage_exons_Kowa")


write.xlsx2(MZFiltered1percentmerge_cov2,"Genes_MiasZapa_reflyrata_arenosa.xlsx",sheetName="Classical_genome_scans",col.names=TRUE,row.names=FALSE)
write.xlsx2(MZ_HighSNP,"Genes_MiasZapa_reflyrata_arenosa.xlsx",sheetName="High_effect_SNPs",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(MZ_INDEL,"Genes_MiasZapa_reflyrata_arenosa.xlsx",sheetName="Indels",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(Mias_cn,"Genes_MiasZapa_reflyrata_arenosa.xlsx",sheetName="CopyNumberVariations",col.names=TRUE,row.names=FALSE,append=T)

write.xlsx2(KKFiltered1percentmerge_cov2,"Genes_KletKowa_reflyrata_arenosa.xlsx",sheetName="Classical_genome_scans",col.names=TRUE,row.names=FALSE)
write.xlsx2(KK_HighSNP,"Genes_KletKowa_reflyrata_arenosa.xlsx",sheetName="High_effect_SNPs",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(KK_INDEL,"Genes_KletKowa_reflyrata_arenosa.xlsx",sheetName="Indels",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(Klet_cn,"Genes_KletKowa_reflyrata_arenosa.xlsx",sheetName="CopyNumberVariations",col.names=TRUE,row.names=FALSE,append=T)



