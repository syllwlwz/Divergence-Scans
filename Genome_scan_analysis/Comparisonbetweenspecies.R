require(xlsx)

#halleri
KK_classical_halleri<-read.xlsx("Genes_KletKowa_reflyrata.xlsx",1,header=T)
KK_high_halleri<-read.xlsx("Genes_KletKowa_reflyrata.xlsx",2,header=T)
KK_Indels_halleri<-read.xlsx("Genes_KletKowa_reflyrata.xlsx",3,header=T)

MZ_classical_halleri<-read.xlsx("Genes_MiasZako_reflyrata.xlsx",1,header=T)
MZ_high_halleri<-read.xlsx("Genes_MiasZako_reflyrata.xlsx",2,header=T)
MZ_Indels_halleri<-read.xlsx("Genes_MiasZako_reflyrata.xlsx",3,header=T)

#arenosa
KK_classical_arenosa<-read.xlsx("Genes_KletKowa_reflyrata_arenosa.xlsx",1,header=T)
KK_high_arenosa<-read.xlsx("Genes_KletKowa_reflyrata_arenosa.xlsx",2,header=T)
KK_Indels_arenosa<-read.xlsx("Genes_KletKowa_reflyrata_arenosa.xlsx",3,header=T)

MZ_classical_arenosa<-read.xlsx("Genes_MiasZapa_reflyrata_arenosa.xlsx",1,header=T)
MZ_high_arenosa<-read.xlsx("Genes_MiasZapa_reflyrata_arenosa.xlsx",2,header=T)
MZ_Indels_arenosa<-read.xlsx("Genes_MiasZapa_reflyrata_arenosa.xlsx",3,header=T)

dup<-duplicated(KK_classical_halleri$Lyr_Gene)
KK_classical_halleri$Rank_in_outlier_metric<-as.numeric(as.character(KK_classical_halleri$Rank_in_outlier_metric))
duplist<-unique(KK_classical_halleri$Lyr_Gene[dup])
KK_classical_halleri_uniquewindows<-KK_classical_halleri[!KK_classical_halleri$Lyr_Gene%in%duplist,]
for (i in 1:length(duplist))
	{Lyrgene<-KK_classical_halleri[KK_classical_halleri$Lyr_Gene==duplist[i],]
	KK_classical_halleri_uniquewindows<-rbind(KK_classical_halleri_uniquewindows,Lyrgene[Lyrgene$Rank_in_outlier_metric==min(Lyrgene$Rank_in_outlier_metric),])
	}

dup<-duplicated(MZ_classical_halleri$Lyr_Gene)
MZ_classical_halleri$Rank_in_outlier_metric<-as.numeric(as.character(MZ_classical_halleri$Rank_in_outlier_metric))
duplist<-unique(MZ_classical_halleri$Lyr_Gene[dup])
MZ_classical_halleri_uniquewindows<-MZ_classical_halleri[!MZ_classical_halleri$Lyr_Gene%in%duplist,]
for (i in 1:length(duplist))
	{Lyrgene<-MZ_classical_halleri[MZ_classical_halleri$Lyr_Gene==duplist[i],]
	MZ_classical_halleri_uniquewindows<-rbind(MZ_classical_halleri_uniquewindows,Lyrgene[Lyrgene$Rank_in_outlier_metric==min(Lyrgene$Rank_in_outlier_metric),])
	}
dup<-duplicated(KK_classical_arenosa$Lyr_Gene)
KK_classical_arenosa$Rank_in_outlier_metric<-as.numeric(as.character(KK_classical_arenosa$Rank_in_outlier_metric))
duplist<-unique(KK_classical_arenosa$Lyr_Gene[dup])
KK_classical_arenosa_uniquewindows<-KK_classical_arenosa[!KK_classical_arenosa$Lyr_Gene%in%duplist,]
for (i in 1:length(duplist))
	{Lyrgene<-KK_classical_arenosa[KK_classical_arenosa$Lyr_Gene==duplist[i],]
	KK_classical_arenosa_uniquewindows<-rbind(KK_classical_arenosa_uniquewindows,Lyrgene[Lyrgene$Rank_in_outlier_metric==min(Lyrgene$Rank_in_outlier_metric),])
	}

dup<-duplicated(MZ_classical_arenosa$Lyr_Gene)
MZ_classical_arenosa$Rank_in_outlier_metric<-as.numeric(as.character(MZ_classical_arenosa$Rank_in_outlier_metric))
duplist<-unique(MZ_classical_arenosa$Lyr_Gene[dup])
MZ_classical_arenosa_uniquewindows<-MZ_classical_arenosa[!MZ_classical_arenosa$Lyr_Gene%in%duplist,]
for (i in 1:length(duplist))
	{Lyrgene<-MZ_classical_arenosa[MZ_classical_arenosa$Lyr_Gene==duplist[i],]
	MZ_classical_arenosa_uniquewindows<-rbind(MZ_classical_arenosa_uniquewindows,Lyrgene[Lyrgene$Rank_in_outlier_metric==min(Lyrgene$Rank_in_outlier_metric),])
	}

KK_classical_halleri_overlap<-KK_classical_halleri_uniquewindows[,c(1:7,73)]
KK_high_halleri_overlap<-KK_high_halleri[,c(7,1:6)]
KK_Indels_halleri_overlap<-KK_Indels_halleri[,c(7,1:6)]

MZ_classical_halleri_overlap<-MZ_classical_halleri_uniquewindows[,c(1:7,73)]
MZ_high_halleri_overlap<-MZ_high_halleri[,c(7,1:6)]
MZ_Indels_halleri_overlap<-MZ_Indels_halleri[,c(7,1:6)]

KK_classical_arenosa_overlap<-KK_classical_arenosa_uniquewindows[,c(1:7,69)]
KK_high_arenosa_overlap<-KK_high_arenosa[,c(7,1:6)]
KK_Indels_arenosa_overlap<-KK_Indels_arenosa[,c(7,1:6)]

MZ_classical_arenosa_overlap<-MZ_classical_arenosa_uniquewindows[,c(1:7,69)]
MZ_high_arenosa_overlap<-MZ_high_arenosa[,c(7,1:6)]
MZ_Indels_arenosa_overlap<-MZ_Indels_arenosa[,c(7,1:6)]


KK_high_halleri_overlap$Outlier_metric<-rep("High_effect_SNP",nrow(KK_high_halleri_overlap))
KK_Indels_halleri_overlap$Outlier_metric<-rep("Indel",nrow(KK_Indels_halleri_overlap))

MZ_high_halleri_overlap$Outlier_metric<-rep("High_effect_SNP",nrow(MZ_high_halleri_overlap))
MZ_Indels_halleri_overlap$Outlier_metric<-rep("Indel",nrow(MZ_Indels_halleri_overlap))

KK_high_arenosa_overlap$Outlier_metric<-rep("High_effect_SNP",nrow(KK_high_arenosa_overlap))
KK_Indels_arenosa_overlap$Outlier_metric<-rep("Indel",nrow(KK_Indels_arenosa_overlap))

MZ_high_arenosa_overlap$Outlier_metric<-rep("High_effect_SNP",nrow(MZ_high_arenosa_overlap))
MZ_Indels_arenosa_overlap$Outlier_metric<-rep("Indel",nrow(MZ_Indels_arenosa_overlap))

KK_halleri_overlap<-rbind(KK_classical_halleri_overlap,KK_high_halleri_overlap,KK_Indels_halleri_overlap)
MZ_halleri_overlap<-rbind(MZ_classical_halleri_overlap,MZ_high_halleri_overlap,MZ_Indels_halleri_overlap)
KK_arenosa_overlap<-rbind(KK_classical_arenosa_overlap,KK_high_arenosa_overlap,KK_Indels_arenosa_overlap)
MZ_arenosa_overlap<-rbind(MZ_classical_arenosa_overlap,MZ_high_arenosa_overlap,MZ_Indels_arenosa_overlap)

KK_halleri_arenosa<-merge(KK_halleri_overlap,KK_arenosa_overlap[,c(1,8)],by="Lyr_Gene")
KK_halleri_arenosa<-KK_halleri_arenosa[!duplicated(KK_halleri_arenosa[,c(1:2,7:9)]),]
names(KK_halleri_arenosa)[8:9]<-c("Outlier_metric_halleri","Outlier_metric_arenosa")

write.table(KK_halleri_arenosa,"KK_halleri_arenosa_overlap_lyrgene.table",sep="\t",row.names=F)

MZ_halleri_arenosa<-merge(MZ_halleri_overlap,MZ_arenosa_overlap[,c(1,8)],by="Lyr_Gene")
MZ_halleri_arenosa<-MZ_halleri_arenosa[!duplicated(MZ_halleri_arenosa[,c(1:2,7:9)]),]
names(MZ_halleri_arenosa)[8:9]<-c("Outlier_metric_halleri","Outlier_metric_arenosa")

write.table(MZ_halleri_arenosa,"MZ_halleri_arenosa_overlap_lyrgene.table",sep="\t",row.names=F)

all<-merge(KK_halleri_arenosa,MZ_halleri_arenosa,by="Lyr_Gene")

KK_MZ_lyr<-read.table("KK_MZ_overlap_lyrgene.table",sep="\t",header=T)

test<-Reduce(function(...) merge(..., by="Lyr_Gene",all.x=TRUE), list(KK_MZ_lyr,KK_high_arenosa,KK_Indels_arenosa,KK_classical_arenosa_uniquewindows,MZ_high_arenosa,MZ_Indels_arenosa,MZ_classical_arenosa_uniquewindows))

#add Araport info
Araport<-read.table("AL_AT_GI_OG_Araport11.txt",sep="\t",header=T,fill=T,quote="")

KK_MZ_arenosa_Araport<-merge(Araport,test,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)

write.table(KK_MZ_arenosa_Araport,"KK_MZ_arenosa_overlap_lyratagenes_filtered0_1perc_Araport.table",sep="\t",row.names=F)

KK_MZ_lyr_halleri<-read.table("KK_MZ_overlap_unique.table.txt",sep="\t",header=T)

KK_MZ_lyr_halleri_merge<-Reduce(function(...) merge(..., by="Lyr_Gene",all.x=TRUE), list(KK_MZ_lyr_halleri,KK_high_halleri,KK_Indels_halleri,KK_classical_halleri_uniquewindows,MZ_high_halleri,MZ_Indels_halleri,MZ_classical_halleri_uniquewindows))

#add Araport info
Araport<-read.table("AL_AT_GI_OG_Araport11.txt",sep="\t",header=T,fill=T,quote="")

KK_MZ_halleri_Araport<-merge(Araport,KK_MZ_lyr_halleri_merge,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)

write.table(KK_MZ_halleri_Araport,"KK_MZ_halleri_overlap_lyratagenes_filtered0_1perc_Araport.table",sep="\t",row.names=F)

MZ_halleri_arenosa_merge<-Reduce(function(...) merge(..., by="Lyr_Gene",all.x=TRUE), list(MZ_halleri_arenosa,MZ_high_halleri,MZ_Indels_halleri,MZ_classical_halleri_uniquewindows,MZ_high_arenosa,MZ_Indels_arenosa,MZ_classical_arenosa_uniquewindows))

#add Araport info
Araport<-read.table("AL_AT_GI_OG_Araport11.txt",sep="\t",header=T,fill=T,quote="")

MZ_halleri_arenosa_Araport<-merge(Araport,MZ_halleri_arenosa_merge,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)

write.table(MZ_halleri_arenosa_Araport,"MZ_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araport.table",sep="\t",row.names=F)

KK_halleri_arenosa_merge<-Reduce(function(...) merge(..., by="Lyr_Gene",all.x=TRUE), list(KK_halleri_arenosa,KK_high_halleri,KK_Indels_halleri,KK_classical_halleri_uniquewindows,KK_high_arenosa,KK_Indels_arenosa,KK_classical_arenosa_uniquewindows))

#add Araport info
Araport<-read.table("AL_AT_GI_OG_Araport11.txt",sep="\t",header=T,fill=T,quote="")

KK_halleri_arenosa_Araport<-merge(Araport,KK_halleri_arenosa_merge,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)

write.table(KK_halleri_arenosa_Araport,"KK_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araport2.table",sep="\t",row.names=F)

template<-read.table("MZ_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araport.table",sep="\t",header=T,check.names=F)

HM_list<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",2,header=T)
HM_list$AGI_code<-toupper(as.character(HM_list$AGI_code))
template_HM<-merge(template,HM_list,by.x="AT_ID",by.y="AGI_code",all.x=T,no.dups=F)
template_HM$AT_ID<-as.character(template_HM$AT_ID)
template_HM$AL_ID<-as.character(template_HM$AL_ID)
geneinfo<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",header=F)
names(geneinfo)<-c("Scaffold","Software","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
geneinfo<-geneinfo[geneinfo$Type=="gene",]
ID<-strsplit(as.character(geneinfo$ID),"=",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T))
geneinfo$ID<-as.character(ID2$X3)
geneinfo$ID<-as.character(geneinfo$ID)
ID<-strsplit(as.character(geneinfo$ID),";",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T))
geneinfo$ID<-as.character(ID2$X1)
template_HM_geneinfo<-merge(template_HM,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
MZ_halleri_arenosa_final<-template_HM_geneinfo[,c(1:11,200:207,187:199,12:186)]
nametemplate<-names(MZ_halleri_arenosa_final)
write.table(MZ_halleri_arenosa_final,"MZ_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded.table",row.names=F,sep="\t")

Mias_coverage_halleri<-read.table("Mias_coverage_genes_exons.table",sep="\t",header=T)
Zako_coverage_halleri<-read.table("Zako_coverage_genes_exons.table",sep="\t",header=T)
Klet_coverage_halleri<-read.table("Klet_coverage_genes_exons.table",sep="\t",header=T)
Kowa_coverage_halleri<-read.table("Kowa_coverage_genes_exons.table",sep="\t",header=T)

Mias_coverage_arenosa<-read.table("Mias_arenosa_genes_exons.table",sep="\t",header=T)
Zapa_coverage_arenosa<-read.table("Zapa_coverage_genes_exons.table",sep="\t",header=T)
Klet_coverage_arenosa<-read.table("Klet_coverage_genes_exons.table",sep="\t",header=T)
Kowa_coverage_arenosa<-read.table("Kowa_coverage_genes_exons.table",sep="\t",header=T)

##############################################################################
#merge all populations#
##############################################################################
MZ_halleri_arenosa_final<-read.table("MZ_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded.table",header=T,sep="\t",check.names=F)

Coverage_halleri<-cbind(Mias_coverage_halleri,Zako_coverage_halleri[,10:13],Klet_coverage_halleri[,10:13],Kowa_coverage_halleri[,10:13])
names(Coverage_halleri)[10:25]<-c("Mias_mean_coverage_gene_halleri","Mias_median_coverage_gene_halleri","Mias_mean_coverage_exons_halleri","Mias_median_coverage_exons_halleri","Zako_mean_coverage_gene_halleri","Zako_median_coverage_gene_halleri","Zako_mean_coverage_exons_halleri","Zako_median_coverage_exons_halleri",
"Klet_mean_coverage_gene_halleri","Klet_median_coverage_gene_halleri","Klet_mean_coverage_exons_halleri","Klet_median_coverage_exons_halleri","Kowa_mean_coverage_gene_halleri","Kowa_median_coverage_gene_halleri","Kowa_mean_coverage_exons_halleri","Kowa_median_coverage_exons_halleri")
Numberofindividuals<-c(8,8,8,8,9,9,9,9,10,10,10,10,7,7,7,7)
Coverage_halleri[,10:25]<-Coverage_halleri[,10:25]/Numberofindividuals

Coverage_arenosa<-cbind(Mias_coverage_arenosa,Zapa_coverage_arenosa[,10:13],Klet_coverage_arenosa[,10:13],Kowa_coverage_arenosa[,10:13])
names(Coverage_arenosa)[10:25]<-c("Mias_mean_coverage_gene_arenosa","Mias_median_coverage_gene_arenosa","Mias_mean_coverage_exons_arenosa","Mias_median_coverage_exons_arenosa","Zapa_mean_coverage_gene_arenosa","Zapa_median_coverage_gene_arenosa","Zapa_mean_coverage_exons_arenosa","Zapa_median_coverage_exons_arenosa",
"Klet_mean_coverage_gene_arenosa","Klet_median_coverage_gene_arenosa","Klet_mean_coverage_exons_arenosa","Klet_median_coverage_exons_arenosa","Kowa_mean_coverage_gene_arenosa","Kowa_median_coverage_gene_arenosa","Kowa_mean_coverage_exons_arenosa","Kowa_median_coverage_exons_arenosa")
Numberofindividuals<-c(28,28,28,28,32,32,32,32,36,36,36,36,32,32,32,32)
Coverage_arenosa[,10:25]<-Coverage_arenosa[,10:25]/Numberofindividuals


Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)
Ahal_ortholist=read.delim("HalleriGeneOrthogroup.txt", header=TRUE)

TH_ortho_numbers<-table(Thal_ortholist$OrthoGroup)
TH_ortho_number<-data.frame(TH_ortho_numbers)
names(TH_ortho_number)<-c("OrthoGroup","Frequency")
Ah_ortho_numbers<-table(Ahal_ortholist$OrthoGroup)
Ah_ortho_number<-data.frame(Ah_ortho_numbers)
names(Ah_ortho_number)<-c("OrthoGroup","Frequency")
Ly_ortho_numbers<-table(Lyr_ortholist$OrthoGroup)
Ly_ortho_number<-data.frame(Ly_ortho_numbers)
names(Ly_ortho_number)<-c("OrthoGroup","Frequency")

MZ_halleri_arenosa_cov<-merge(MZ_halleri_arenosa_final[,1:4],TH_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)
names(MZ_halleri_arenosa_cov)[5]<-"Thaliana_paralogue_number"
MZ_halleri_arenosa_cov2<-merge(MZ_halleri_arenosa_cov,Ah_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)
names(MZ_halleri_arenosa_cov2)[6]<-"Halleri_paralogue_number"
MZ_halleri_arenosa_cov3<-merge(MZ_halleri_arenosa_cov2,Ly_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)
names(MZ_halleri_arenosa_cov3)[7]<-"Lyrata_paralogue_number"

MZ_halleri_arenosa_cov_halleri<-merge(MZ_halleri_arenosa_final[,1:4],Coverage_halleri,by.x="AL_ID",by.y="ID",all.x=T)
MZ_halleri_arenosa_cov_arenosa<-merge(MZ_halleri_arenosa_final[,1:4],Coverage_arenosa,by.x="AL_ID",by.y="ID",all.x=T)

MZ_halleri_arenosa_final_wcov<-append(MZ_halleri_arenosa_final,MZ_halleri_arenosa_cov_halleri[,13:28],after=34)
MZ_halleri_arenosa_final_wcov2<-append(MZ_halleri_arenosa_final_wcov,MZ_halleri_arenosa_cov_arenosa[,13:28],after=50)
MZ_halleri_arenosa_final_wcov3<-append(MZ_halleri_arenosa_final_wcov2,MZ_halleri_arenosa_cov3[,5:7],after=66)

write.table(MZ_halleri_arenosa_final_wcov3,"MZ_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded.table",sep="\t",row.names=F)

KK_halleri_arenosa_final<-read.table("KK_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araport2.table",header=T,sep="\t",check.names=F)
HM_list<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",2,header=T)
HM_list$AGI_code<-toupper(as.character(HM_list$AGI_code))
KK_halleri_arenosa_final_HM<-merge(KK_halleri_arenosa_final,HM_list,by.x="AT_ID",by.y="AGI_code",all.x=T,no.dups=F)
KK_halleri_arenosa_final_HM$AT_ID<-as.character(KK_halleri_arenosa_final_HM$AT_ID)
KK_halleri_arenosa_final_HM$AL_ID<-as.character(KK_halleri_arenosa_final_HM$AL_ID)
geneinfo<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",header=F)
names(geneinfo)<-c("Scaffold","Software","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
geneinfo<-geneinfo[geneinfo$Type=="gene",]
ID<-strsplit(as.character(geneinfo$ID),"=",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T))
geneinfo$ID<-as.character(ID2$X3)
geneinfo$ID<-as.character(geneinfo$ID)
ID<-strsplit(as.character(geneinfo$ID),";",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T))
geneinfo$ID<-as.character(ID2$X1)
KK_halleri_arenosa_final_HM_geneinfo<-merge(KK_halleri_arenosa_final_HM,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_halleri_arenosa_final2<-KK_halleri_arenosa_final_HM_geneinfo[,c(1:11,340:347,327:339,18:19,319:326,26:319)]

TH_ortho_numbers<-table(Thal_ortholist$OrthoGroup)
TH_ortho_number<-data.frame(TH_ortho_numbers)
names(TH_ortho_number)<-c("OrthoGroup","Frequency")
Ah_ortho_numbers<-table(Ahal_ortholist$OrthoGroup)
Ah_ortho_number<-data.frame(Ah_ortho_numbers)
names(Ah_ortho_number)<-c("OrthoGroup","Frequency")
Ly_ortho_numbers<-table(Lyr_ortholist$OrthoGroup)
Ly_ortho_number<-data.frame(Ly_ortho_numbers)
names(Ly_ortho_number)<-c("OrthoGroup","Frequency")

KK_halleri_arenosa_cov<-merge(KK_halleri_arenosa_final2[,1:4],TH_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)
names(KK_halleri_arenosa_cov)[5]<-"Thaliana_paralogue_number"
KK_halleri_arenosa_cov2<-merge(KK_halleri_arenosa_cov,Ah_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)
names(KK_halleri_arenosa_cov2)[6]<-"Halleri_paralogue_number"
KK_halleri_arenosa_cov3<-merge(KK_halleri_arenosa_cov2,Ly_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)
names(KK_halleri_arenosa_cov3)[7]<-"Lyrata_paralogue_number"

KK_halleri_arenosa_cov_halleri<-merge(KK_halleri_arenosa_final2[,1:4],Coverage_halleri,by.x="AL_ID",by.y="ID",all.x=T)
KK_halleri_arenosa_cov_arenosa<-merge(KK_halleri_arenosa_final2[,1:4],Coverage_arenosa,by.x="AL_ID",by.y="ID",all.x=T)

KK_halleri_arenosa_final_wcov<-append(KK_halleri_arenosa_final2,KK_halleri_arenosa_cov_halleri[,13:28],after=34)
KK_halleri_arenosa_final_wcov2<-append(KK_halleri_arenosa_final_wcov,KK_halleri_arenosa_cov_arenosa[,13:28],after=50)
KK_halleri_arenosa_final_wcov3<-append(KK_halleri_arenosa_final_wcov2,KK_halleri_arenosa_cov3[,5:7],after=66)
KK_halleri_arenosa_final_wcov4<-as.data.frame(KK_halleri_arenosa_final_wcov3)
KK_halleri_arenosa_final_wcov5<-KK_halleri_arenosa_final_wcov4[,-c(70:77)]
KK_halleri_arenosa_final_wcov6<-KK_halleri_arenosa_final_wcov5[,-c(77:94)]
KK_halleri_arenosa_final_wcov7<-KK_halleri_arenosa_final_wcov6[,-c(79:84)]
KK_halleri_arenosa_final_wcov8<-KK_halleri_arenosa_final_wcov7[,-c(88:100)]
KK_halleri_arenosa_final_wcov9<-KK_halleri_arenosa_final_wcov8[,-c(100:105)]
KK_halleri_arenosa_final_wcov10<-KK_halleri_arenosa_final_wcov9[,-c(122:126)]
KK_halleri_arenosa_final_wcov11<-KK_halleri_arenosa_final_wcov10[,-c(147:159)]
KK_halleri_arenosa_final_wcov12<-KK_halleri_arenosa_final_wcov11[,-c(157:162)]
KK_halleri_arenosa_final_wcov13<-KK_halleri_arenosa_final_wcov12[,-c(164:181)]
KK_halleri_arenosa_final_wcov14<-KK_halleri_arenosa_final_wcov13[,-c(169:174)]
KK_halleri_arenosa_final_wcov15<-KK_halleri_arenosa_final_wcov14[,-c(178:190)]
KK_halleri_arenosa_final_wcov16<-KK_halleri_arenosa_final_wcov15[,-c(190:195)]
KK_halleri_arenosa_final_wcov17<-KK_halleri_arenosa_final_wcov16[,-c(208:212)]
KK_halleri_arenosa_final_wcov18<-KK_halleri_arenosa_final_wcov17[,-c(233:245)]

write.table(names(MZ_halleri_arenosa_final_wcov3),"Nametemplate.table",sep="\t",row.names=F)
Nametemplate<-read.table("Nametemplate.table",sep="\t",header=F)
names(KK_halleri_arenosa_final_wcov18)<-Nametemplate$V1[1:235]
write.table(KK_halleri_arenosa_final_wcov18,"KK_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded.table",sep="\t",row.names=F)



KK_MZ_halleri<-read.table("KK_MZ_halleri_overlap_lyratagenes_filtered0_1perc_Araport.table.txt",header=T,sep="\t",check.names=F)
HM_list<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",2,header=T)
HM_list$AGI_code<-toupper(as.character(HM_list$AGI_code))
KK_MZ_halleri_HM<-merge(KK_MZ_halleri,HM_list,by.x="AT_ID",by.y="AGI_code",all.x=T,no.dups=F)
KK_MZ_halleri_HM$AT_ID<-as.character(KK_MZ_halleri_HM$AT_ID)
KK_MZ_halleri_HM$AL_ID<-as.character(KK_MZ_halleri_HM$AL_ID)
geneinfo<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",header=F)
names(geneinfo)<-c("Scaffold","Software","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
geneinfo<-geneinfo[geneinfo$Type=="gene",]
ID<-strsplit(as.character(geneinfo$ID),"=",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T))
geneinfo$ID<-as.character(ID2$X3)
geneinfo$ID<-as.character(geneinfo$ID)
ID<-strsplit(as.character(geneinfo$ID),";",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T))
geneinfo$ID<-as.character(ID2$X1)
KK_MZ_halleri_HM_geneinfo<-merge(KK_MZ_halleri_HM,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_MZ_halleri2<-KK_MZ_halleri_HM_geneinfo[,c(1:11,526:533,513:525,18,25,504:511,26:503)]

TH_ortho_numbers<-table(Thal_ortholist$OrthoGroup)
TH_ortho_number<-data.frame(TH_ortho_numbers)
names(TH_ortho_number)<-c("OrthoGroup","Frequency")
Ah_ortho_numbers<-table(Ahal_ortholist$OrthoGroup)
Ah_ortho_number<-data.frame(Ah_ortho_numbers)
names(Ah_ortho_number)<-c("OrthoGroup","Frequency")
Ly_ortho_numbers<-table(Lyr_ortholist$OrthoGroup)
Ly_ortho_number<-data.frame(Ly_ortho_numbers)
names(Ly_ortho_number)<-c("OrthoGroup","Frequency")

KK_MZ_halleri_cov<-merge(KK_MZ_halleri2[,1:4],TH_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)
names(KK_MZ_halleri_cov)[5]<-"Thaliana_paralogue_number"
KK_MZ_halleri_cov2<-merge(KK_MZ_halleri_cov,Ah_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)
names(KK_MZ_halleri_cov2)[6]<-"Halleri_paralogue_number"
KK_MZ_halleri_cov3<-merge(KK_MZ_halleri_cov2,Ly_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)
names(KK_MZ_halleri_cov3)[7]<-"Lyrata_paralogue_number"

KK_MZ_halleri_cov_halleri<-merge(KK_MZ_halleri2[,1:4],Coverage_halleri,by.x="AL_ID",by.y="ID",all.x=T)

KK_MZ_halleri_wcov<-append(KK_MZ_halleri2,KK_MZ_halleri_cov_halleri[,13:28],after=34)
KK_MZ_halleri_wcov3<-append(KK_MZ_halleri_wcov,KK_MZ_halleri_cov3[,5:7],after=58)
KK_MZ_halleri_wcov4<-as.data.frame(KK_MZ_halleri_wcov3)
KK_MZ_halleri_wcov5<-KK_MZ_halleri_wcov4[,-c(69:81)]
KK_MZ_halleri_wcov6<-KK_MZ_halleri_wcov5[,-c(240:257)]
KK_MZ_halleri_wcov7<-KK_MZ_halleri_wcov6[,-c(227:232)]
KK_MZ_halleri_wcov8<-KK_MZ_halleri_wcov7[,-c(251:263)]
KK_MZ_halleri_wcov9<-KK_MZ_halleri_wcov8[,-c(236:241)]
KK_MZ_halleri_wcov10<-KK_MZ_halleri_wcov9[,-c(257:262)]
KK_MZ_halleri_wcov11<-KK_MZ_halleri_wcov10[,-c(279:283)]
KK_MZ_halleri_wcov12<-KK_MZ_halleri_wcov11[,-c(304:316)]
KK_MZ_halleri_wcov13<-KK_MZ_halleri_wcov12[,-c(314:319)]
KK_MZ_halleri_wcov14<-KK_MZ_halleri_wcov13[,-c(321:338)]
KK_MZ_halleri_wcov15<-KK_MZ_halleri_wcov14[,-c(323:328)]
KK_MZ_halleri_wcov16<-KK_MZ_halleri_wcov15[,-c(332:344)]
KK_MZ_halleri_wcov17<-KK_MZ_halleri_wcov16[,-c(344:349)]
KK_MZ_halleri_wcov18<-KK_MZ_halleri_wcov17[,-c(366:370)]
KK_MZ_halleri_wcov19<-KK_MZ_halleri_wcov18[,-c(391:403)]

write.table(KK_MZ_halleri_wcov19,"KK_MZ_halleri_overlap_lyratagenes_filtered0_1perc_Araportadded.table",sep="\t",row.names=F)

KK_MZ_halleri_wcov19<-read.table("KK_MZ_halleri_overlap_lyratagenes_filtered0_1perc_Araportadded.csv",sep=";",header=T)




KK_MZ_arenosa<-read.table("KK_MZ_arenosa_overlap_lyratagenes_filtered0_1perc_Araport.table.txt",header=T,sep="\t",check.names=F)

KK_MZ_arenosa_HM_geneinfo<-merge(KK_MZ_arenosa,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_MZ_arenosa2<-KK_MZ_arenosa_HM_geneinfo[,c(1:11,212:219,20:33,18:19,203:210,34:202)]

TH_ortho_numbers<-table(Thal_ortholist$OrthoGroup)
TH_ortho_number<-data.frame(TH_ortho_numbers)
names(TH_ortho_number)<-c("OrthoGroup","Frequency")
Ah_ortho_numbers<-table(Ahal_ortholist$OrthoGroup)
Ah_ortho_number<-data.frame(Ah_ortho_numbers)
names(Ah_ortho_number)<-c("OrthoGroup","Frequency")
Ly_ortho_numbers<-table(Lyr_ortholist$OrthoGroup)
Ly_ortho_number<-data.frame(Ly_ortho_numbers)
names(Ly_ortho_number)<-c("OrthoGroup","Frequency")

KK_MZ_arenosa_cov<-merge(KK_MZ_arenosa2[,1:4],TH_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)
names(KK_MZ_arenosa_cov)[5]<-"Thaliana_paralogue_number"
KK_MZ_arenosa_cov2<-merge(KK_MZ_arenosa_cov,Ah_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)
names(KK_MZ_arenosa_cov2)[6]<-"Halleri_paralogue_number"
KK_MZ_arenosa_cov3<-merge(KK_MZ_arenosa_cov2,Ly_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)
names(KK_MZ_arenosa_cov3)[7]<-"Lyrata_paralogue_number"

KK_MZ_arenosa_cov_arenosa<-merge(KK_MZ_arenosa2[,1:4],Coverage_arenosa,by.x="AL_ID",by.y="ID",all.x=T)

KK_MZ_arenosa_wcov<-append(KK_MZ_arenosa2,KK_MZ_arenosa_cov_arenosa[,13:28],after=35)
KK_MZ_arenosa_wcov3<-append(KK_MZ_arenosa_wcov,KK_MZ_arenosa_cov3[,5:7],after=59)
KK_MZ_arenosa_wcov4<-as.data.frame(KK_MZ_arenosa_wcov3)
KK_MZ_arenosa_wcov5<-KK_MZ_arenosa_wcov4[,-c(21)]
KK_MZ_arenosa_wcov6<-KK_MZ_arenosa_wcov5[,-c(69:73)]
KK_MZ_arenosa_wcov7<-KK_MZ_arenosa_wcov6[,-c(71:73)]
KK_MZ_arenosa_wcov8<-KK_MZ_arenosa_wcov7[,-c(137:144)]
KK_MZ_arenosa_wcov9<-KK_MZ_arenosa_wcov8[,-c(146:148)]

template<-names(KK_MZ_halleri_wcov19)[-c(110:114)]
template2<-c(template[1:128],"Rank_in_outlier_metric",template[129:214])
template3<-template[-c(185:188)]
template4<-c(template3[1:203],"Rank_in_outlier_metric",template3[204:210])

names(KK_MZ_arenosa_wcov9)<-template4

write.table(KK_MZ_arenosa_wcov9,"KK_MZ_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded.table",sep="\t",row.names=F)


#new Araport transformation of downloaded file because of error in old and new combination of tables with Araport annotation
Araport_description<-read.table("Araport11_GFF3_genes_transposons.201606.gff",sep="\t",header=F,fill=T,quote="")
Araport_description<-Araport_description[Araport_description$V3=="gene",]
write.table(Araport_description,"Araport11_GFF3_genes_transposons.201606_V2.gff",sep="\t",row.names=F)
ID<-strsplit(as.character(Araport_description$V9),";",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T,ncol=10))
ID3<-lapply(ID2,2,function(x) as.character(x))
ID4<-data.frame(ID3,stringsAsFactors=FALSE)

require(qdap)
locus_type<-gsub("NA","",paste2(matrix(ID4[apply(ID4[1:nrow(ID4),],1:2,function(x) grep("locus_type=",x))==1],nrow=nrow(ID4),ncol=10),sep="",handle.na=F))
locus_type2<-unlist(lapply(strsplit(as.character(locus_type),"=",fixed=T),"[",2))
symbol<-gsub("NA","",paste2(matrix(ID4[apply(ID4[1:nrow(ID4),],1:2,function(x) grep("symbol=",x))==1],nrow=nrow(ID4),ncol=10),sep="",handle.na=F))
symbol2<-unlist(lapply(strsplit(as.character(symbol),"=",fixed=T),"[",2))
AT_ID<-gsub("NA","",paste2(matrix(ID4[apply(ID4[1:nrow(ID4),],1:2,function(x) grep("ID=",x))==1],nrow=nrow(ID4),ncol=10),sep="",handle.na=F))
AT_ID2<-unlist(lapply(strsplit(as.character(AT_ID),"=",fixed=T),"[",2))
note<-gsub("NA","",paste2(matrix(ID4[apply(ID4[1:nrow(ID4),],1:2,function(x) grep("Note=",x))==1],nrow=nrow(ID4),ncol=10),sep="",handle.na=F))
note2<-unlist(lapply(strsplit(as.character(note),"=",fixed=T),"[",2))
full_name<-gsub("NA","",paste2(matrix(ID4[apply(ID4[1:nrow(ID4),],1:2,function(x) grep("full_name=",x))==1],nrow=nrow(ID4),ncol=10),sep="",handle.na=F))
full_name2<-unlist(lapply(strsplit(as.character(full_name),"=",fixed=T),"[",2))
curator_summary<-gsub("NA","",paste2(matrix(ID4[apply(ID4[1:nrow(ID4),],1:2,function(x) grep("curator_summary=",x))==1],nrow=nrow(ID4),ncol=10),sep="",handle.na=F))
curator_summary2<-unlist(lapply(strsplit(as.character(curator_summary),"=",fixed=T),"[",2))
description<-gsub("NA","",paste2(matrix(ID4[apply(ID4[1:nrow(ID4),],1:2,function(x) grep("description=",x))==1],nrow=nrow(ID4),ncol=10),sep="",handle.na=F))
description2<-unlist(lapply(strsplit(as.character(description),"=",fixed=T),"[",2))
computational_description<-gsub("NA","",paste2(matrix(ID4[apply(ID4[1:nrow(ID4),],1:2,function(x) grep("computational_description=",x))==1],nrow=nrow(ID4),ncol=10),sep="",handle.na=F))
computational_description2<-unlist(lapply(strsplit(as.character(computational_description),"=",fixed=T),"[",2))

Araport_df<-data.frame(AT_ID2,symbol2,locus_type2,note2,full_name2,description2,curator_summary2,computational_description2)
names(Araport_df)<-c("AT_ID","symbol","locus_type","note","full_name","description","curator_summary","computational_description")

KK_MZ_arenosa<-read.table("KK_MZ_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded.csv",sep=";",header=T)
KK_MZ_halleri<-read.table("KK_MZ_halleri_overlap_lyratagenes_filtered0_1perc_Araportadded.csv",sep=";",header=T)
KK_halleri_arenosa<-read.table("KK_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded.csv",sep=";",header=T)
MZ_halleri_arenosa<-read.table("MZ_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded.csv",sep=";",header=T)

KK_MZ_arenosa$OrthoGroup<-as.character(KK_MZ_arenosa$OrthoGroup)
KK_MZ_halleri$OrthoGroup<-as.character(KK_MZ_halleri$OrthoGroup)
KK_halleri_arenosa$OrthoGroup<-as.character(KK_halleri_arenosa$OrthoGroup)
MZ_halleri_arenosa$OrthoGroup<-as.character(MZ_halleri_arenosa$OrthoGroup)
Thal_ortholist$OrthoGroup<-as.character(Thal_ortholist$OrthoGroup)

KK_MZ_arenosa_thal<-Thal_ortholist[match(KK_MZ_arenosa$OrthoGroup,Thal_ortholist$OrthoGroup),]
KK_MZ_halleri_thal<-Thal_ortholist[match(KK_MZ_halleri$OrthoGroup,Thal_ortholist$OrthoGroup),]
KK_halleri_arenosa_thal<-Thal_ortholist[match(KK_halleri_arenosa$OrthoGroup,Thal_ortholist$OrthoGroup),]
MZ_halleri_arenosa_thal<-Thal_ortholist[match(MZ_halleri_arenosa$OrthoGroup,Thal_ortholist$OrthoGroup),]

Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)

KK_MZ_arenosa_araport<-Araport_df[match(KK_MZ_arenosa_thal$Gene,Araport_df$AT_ID),]
KK_MZ_halleri_araport<-Araport_df[match(KK_MZ_halleri_thal$Gene,Araport_df$AT_ID),]
KK_halleri_arenosa_araport<-Araport_df[match(KK_halleri_arenosa_thal$Gene,Araport_df$AT_ID),]
MZ_halleri_arenosa_araport<-Araport_df[match(MZ_halleri_arenosa_thal$Gene,Araport_df$AT_ID),]

KK_MZ_arenosa_HM<-HM_list[match(KK_MZ_arenosa_thal$Gene,HM_list$AGI_code),]
KK_MZ_halleri_HM<-HM_list[match(KK_MZ_halleri_thal$Gene,HM_list$AGI_code),]
KK_halleri_arenosa_HM<-HM_list[match(KK_halleri_arenosa_thal$Gene,HM_list$AGI_code),]
MZ_halleri_arenosa_HM<-HM_list[match(MZ_halleri_arenosa_thal$Gene,HM_list$AGI_code),]


KK_MZ_arenosa2<-data.frame(KK_MZ_arenosa_araport,KK_MZ_arenosa_HM)
KK_MZ_halleri2<-data.frame(KK_MZ_halleri_araport,KK_MZ_halleri_HM)
KK_halleri_arenosa2<-data.frame(KK_halleri_arenosa_araport,KK_halleri_arenosa_HM)
MZ_halleri_arenosa2<-data.frame(MZ_halleri_arenosa_araport,MZ_halleri_arenosa_HM)

write.table(KK_MZ_arenosa2,"KK_MZ_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded_corrected.table",sep="\t",row.names=F)
write.table(KK_MZ_halleri2,"KK_MZ_halleri_overlap_lyratagenes_filtered0_1perc_Araportadded_corrected.table",sep="\t",row.names=F)
write.table(KK_halleri_arenosa2,"KK_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded_corrected.table",sep="\t",row.names=F)
write.table(MZ_halleri_arenosa2,"MZ_halleri_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded_corrected.table",sep="\t",row.names=F)


Araport_df$AT_ID<-as.character(Araport_df$AT_ID)
Thal_ortholist$Gene<-as.character(Thal_ortholist$Gene)
AL_AT_Araport_df1<-merge(Araport_df,Thal_ortholist,by.x="AT_ID",by.y="Gene",all.x=T)
AL_AT_Araport_df1$OrthoGroup<-as.character(AL_AT_Araport_df1$OrthoGroup)
Lyr_ortholist$OrthoGroup<-as.character(Lyr_ortholist$OrthoGroup)
AL_AT_Araport_df<-merge(AL_AT_Araport_df1,Lyr_ortholist,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)

AL_AT_Araport_df2<-AL_AT_Araport_df[,c(10,2,1,3:9)]
names(AL_AT_Araport_df2)[1]<-"AL_ID"
AL_AT_Araport_df2$AL_ID<-as.character(AL_AT_Araport_df2$AL_ID)


write.table(AL_AT_Araport_df2,"Araport_AL_AT.table",sep="\t",row.names=F)


Araport<-read.table("Araport_AL_AT.table",sep="\t",header=T)
HM_list<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",2,header=T)
HM_list$AGI_code<-toupper(as.character(HM_list$AGI_code))

Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)
Ahal_ortholist=read.delim("HalleriGeneOrthogroup.txt", header=TRUE)

TH_ortho_numbers<-table(Thal_ortholist$OrthoGroup)
TH_ortho_number<-data.frame(TH_ortho_numbers)
names(TH_ortho_number)<-c("OrthoGroup","Frequency")
Ah_ortho_numbers<-table(Ahal_ortholist$OrthoGroup)
Ah_ortho_number<-data.frame(Ah_ortho_numbers)
names(Ah_ortho_number)<-c("OrthoGroup","Frequency")
Ly_ortho_numbers<-table(Lyr_ortholist$OrthoGroup)
Ly_ortho_number<-data.frame(Ly_ortho_numbers)
names(Ly_ortho_number)<-c("OrthoGroup","Frequency")

geneinfo<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",header=F)
names(geneinfo)<-c("Scaffold","Software","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
geneinfo<-geneinfo[geneinfo$Type=="gene",]
ID<-strsplit(as.character(geneinfo$ID),"=",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T))
geneinfo$ID<-as.character(ID2$X3)
geneinfo$ID<-as.character(geneinfo$ID)
ID<-strsplit(as.character(geneinfo$ID),";",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T))
geneinfo$ID<-as.character(ID2$X1)

Mias_coverage_halleri<-read.table("Mias_coverage_genes_exons.table",sep="\t",header=T)
Zako_coverage_halleri<-read.table("Zako_coverage_genes_exons.table",sep="\t",header=T)
Klet_coverage_halleri<-read.table("Klet_coverage_genes_exons.table",sep="\t",header=T)
Kowa_coverage_halleri<-read.table("Kowa_coverage_genes_exons.table",sep="\t",header=T)

Mias_coverage_arenosa<-read.table("Mias_arenosa_genes_exons.table",sep="\t",header=T)
Zapa_coverage_arenosa<-read.table("Zapa_coverage_genes_exons.table",sep="\t",header=T)
Klet_coverage_arenosa<-read.table("Klet_coverage_genes_exons.table",sep="\t",header=T)
Kowa_coverage_arenosa<-read.table("Kowa_coverage_genes_exons.table",sep="\t",header=T)

Coverage_halleri<-cbind(Mias_coverage_halleri,Zako_coverage_halleri[,10:13],Klet_coverage_halleri[,10:13],Kowa_coverage_halleri[,10:13])
names(Coverage_halleri)[10:25]<-c("Mias_mean_coverage_gene_halleri","Mias_median_coverage_gene_halleri","Mias_mean_coverage_exons_halleri","Mias_median_coverage_exons_halleri","Zako_mean_coverage_gene_halleri","Zako_median_coverage_gene_halleri","Zako_mean_coverage_exons_halleri","Zako_median_coverage_exons_halleri",
"Klet_mean_coverage_gene_halleri","Klet_median_coverage_gene_halleri","Klet_mean_coverage_exons_halleri","Klet_median_coverage_exons_halleri","Kowa_mean_coverage_gene_halleri","Kowa_median_coverage_gene_halleri","Kowa_mean_coverage_exons_halleri","Kowa_median_coverage_exons_halleri")
Numberofindividuals<-c(8,8,8,8,9,9,9,9,10,10,10,10,7,7,7,7)
Coverage_halleri[,10:25]<-Coverage_halleri[,10:25]/Numberofindividuals

Coverage_arenosa<-cbind(Mias_coverage_arenosa,Zapa_coverage_arenosa[,10:13],Klet_coverage_arenosa[,10:13],Kowa_coverage_arenosa[,10:13])
names(Coverage_arenosa)[10:25]<-c("Mias_mean_coverage_gene_arenosa","Mias_median_coverage_gene_arenosa","Mias_mean_coverage_exons_arenosa","Mias_median_coverage_exons_arenosa","Zapa_mean_coverage_gene_arenosa","Zapa_median_coverage_gene_arenosa","Zapa_mean_coverage_exons_arenosa","Zapa_median_coverage_exons_arenosa",
"Klet_mean_coverage_gene_arenosa","Klet_median_coverage_gene_arenosa","Klet_mean_coverage_exons_arenosa","Klet_median_coverage_exons_arenosa","Kowa_mean_coverage_gene_arenosa","Kowa_median_coverage_gene_arenosa","Kowa_mean_coverage_exons_arenosa","Kowa_median_coverage_exons_arenosa")
Numberofindividuals<-c(28,28,28,28,32,32,32,32,36,36,36,36,32,32,32,32)
Coverage_arenosa[,10:25]<-Coverage_arenosa[,10:25]/Numberofindividuals



require(xlsx)

MZ_GS_halleri<-read.xlsx2("Genes_MiasZako_reflyrata.xlsx",1,header=T)
MZ_GS_halleri_Araport<-merge(Araport,MZ_GS_halleri,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
MZ_GS_halleri_Araport_geneinfo<-merge(MZ_GS_halleri_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
MZ_GS_halleri_Araport_geneinfo2<-MZ_GS_halleri_Araport_geneinfo[,c(1:11,93:100,18:92)]

MZ_GS_halleri_highSNPs<-read.xlsx2("Genes_MiasZako_reflyrata.xlsx",2,header=T)
MZ_GS_halleri_highSNPs_Araport<-merge(Araport,MZ_GS_halleri_highSNPs,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
MZ_GS_halleri_highSNPs_Araport_geneinfo<-merge(MZ_GS_halleri_highSNPs_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
MZ_GS_halleri_highSNPs_geneinfo2<-MZ_GS_halleri_highSNPs_Araport_geneinfo[,c(1:11,45:52,18:44)]


write.table(MZ_GS_halleri_highSNPs_Araport,"MZ_GS_halleri_Araport_highSNPs.table",sep="\t",row.names=F)

MZ_GS_halleri_Indels<-read.xlsx2("Genes_MiasZako_reflyrata.xlsx",3,header=T)
MZ_GS_halleri_Indels_Araport<-merge(Araport,MZ_GS_halleri_Indels,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
MZ_GS_halleri_Indels_Araport_geneinfo<-merge(MZ_GS_halleri_Indels_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
MZ_GS_halleri_Indels_Araport_geneinfo2<-MZ_GS_halleri_Indels_Araport_geneinfo[,c(1:11,52:59,18:51)]


write.table(MZ_GS_halleri_Indels_Araport,"MZ_GS_halleri_Araport_Indels.table",sep="\t",row.names=F)
write.xlsx(MZ_GS_halleri_Araport,"Genes_MiasZako_halleri_Araportadded.xlsx",sheetName="Classical_genome_scans",col.names=TRUE,row.names=FALSE)
write.xlsx2(MZ_GS_halleri_highSNPs_Araport,"Genes_MiasZako_halleri_Araportadded.xlsx",sheetName="High_effect_SNPs",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(MZ_GS_halleri_Indels_Araport,"Genes_MiasZako_halleri_Araportadded.xlsx",sheetName="Indels",col.names=TRUE,row.names=FALSE,append=T)




KK_GS_halleri<-read.xlsx2("Genes_KletKowa_reflyrata.xlsx",1,header=T)
KK_GS_halleri_Araport<-merge(Araport,KK_GS_halleri,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
KK_GS_halleri_Araport_geneinfo<-merge(KK_GS_halleri_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_GS_halleri_Araport_geneinfo2<-KK_GS_halleri_Araport_geneinfo[,c(1:11,92:99,18:91)]

KK_GS_halleri_highSNPs<-read.xlsx2("Genes_KletKowa_reflyrata.xlsx",2,header=T)
KK_GS_halleri_highSNPs_Araport<-merge(Araport,KK_GS_halleri_highSNPs,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
KK_GS_halleri_highSNPs_Araport_geneinfo<-merge(KK_GS_halleri_highSNPs_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_GS_halleri_highSNPs_geneinfo2<-KK_GS_halleri_highSNPs_Araport_geneinfo[,c(1:11,44:51,18:43)]

KK_GS_halleri_Indels<-read.xlsx2("Genes_KletKowa_reflyrata.xlsx",3,header=T)
KK_GS_halleri_Indels_Araport<-merge(Araport,KK_GS_halleri_Indels,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
KK_GS_halleri_Indels_Araport_geneinfo<-merge(KK_GS_halleri_Indels_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_GS_halleri_Indels_Araport_geneinfo2<-KK_GS_halleri_Indels_Araport_geneinfo[,c(1:11,51:58,18:50)]

write.table(KK_GS_halleri_Indels_Araport,"KK_GS_halleri_Araport_Indels.table",sep="\t",row.names=F)
write.xlsx(KK_GS_halleri_Araport,"Genes_KletKowa_halleri_Araportadded.xlsx",sheetName="Classical_genome_scans",col.names=TRUE,row.names=FALSE)
write.xlsx2(KK_GS_halleri_highSNPs_Araport,"Genes_KletKowa_halleri_Araportadded.xlsx",sheetName="High_effect_SNPs",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(KK_GS_halleri_Indels_Araport,"Genes_KletKowa_halleri_Araportadded.xlsx",sheetName="Indels",col.names=TRUE,row.names=FALSE,append=T)








MZ_GS_arenosa<-read.xlsx2("Genes_MiasZapa_reflyrata_arenosa.xlsx",1,header=T)
MZ_GS_arenosa_Araport<-merge(Araport,MZ_GS_arenosa,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
MZ_GS_arenosa_Araport_geneinfo<-merge(MZ_GS_arenosa_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
MZ_GS_arenosa_Araport_geneinfo2<-MZ_GS_arenosa_Araport_geneinfo[,c(1:10,88:95,17:87)]

MZ_GS_arenosa_highSNPs<-read.xlsx2("Genes_MiasZapa_reflyrata_arenosa.xlsx",2,header=T)
MZ_GS_arenosa_highSNPs_Araport<-merge(Araport,MZ_GS_arenosa_highSNPs,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
MZ_GS_arenosa_highSNPs_Araport_geneinfo<-merge(MZ_GS_arenosa_highSNPs_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
MZ_GS_arenosa_highSNPs_geneinfo2<-MZ_GS_arenosa_highSNPs_Araport_geneinfo[,c(1:10,47:54,17:46)]

MZ_GS_arenosa_Indels<-read.xlsx2("Genes_MiasZapa_reflyrata_arenosa.xlsx",3,header=T)
MZ_GS_arenosa_Indels_Araport<-merge(Araport,MZ_GS_arenosa_Indels,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
MZ_GS_arenosa_Indels_Araport_geneinfo<-merge(MZ_GS_arenosa_Indels_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
MZ_GS_arenosa_Indels_Araport_geneinfo2<-MZ_GS_arenosa_Indels_Araport_geneinfo[,c(1:10,51:58,17:50)]

write.xlsx(MZ_GS_arenosa_Araport,"Genes_MiasZapa_arenosa_Araportadded.xlsx",sheetName="Classical_genome_scans",col.names=TRUE,row.names=FALSE)
write.xlsx2(MZ_GS_arenosa_highSNPs_Araport,"Genes_MiasZapa_arenosa_Araportadded.xlsx",sheetName="High_effect_SNPs",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(MZ_GS_arenosa_Indels_Araport,"Genes_MiasZapa_arenosa_Araportadded.xlsx",sheetName="Indels",col.names=TRUE,row.names=FALSE,append=T)




KK_GS_arenosa<-read.xlsx2("Genes_KletKowa_reflyrata_arenosa.xlsx",1,header=T)
KK_GS_arenosa_Araport<-merge(Araport,KK_GS_arenosa,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
KK_GS_arenosa_Araport_geneinfo<-merge(KK_GS_arenosa_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_GS_arenosa_Araport_geneinfo2<-KK_GS_arenosa_Araport_geneinfo[,c(1:10,88:95,17:87)]

KK_GS_arenosa_highSNPs<-read.xlsx2("Genes_KletKowa_reflyrata_arenosa.xlsx",2,header=T)
KK_GS_arenosa_highSNPs_Araport<-merge(Araport,KK_GS_arenosa_highSNPs,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
KK_GS_arenosa_highSNPs_Araport_geneinfo<-merge(KK_GS_arenosa_highSNPs_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_GS_arenosa_highSNPs_geneinfo2<-KK_GS_arenosa_highSNPs_Araport_geneinfo[,c(1:11,47:54,17:46)]

KK_GS_arenosa_Indels<-read.xlsx2("Genes_KletKowa_reflyrata_arenosa.xlsx",3,header=T)
KK_GS_arenosa_Indels_Araport<-merge(Araport,KK_GS_arenosa_Indels,by.y="Lyr_Gene",by.x="AL_ID",all.y=T)
KK_GS_arenosa_Indels_Araport_geneinfo<-merge(KK_GS_arenosa_Indels_Araport,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_GS_arenosa_Indels_Araport_geneinfo2<-KK_GS_arenosa_Indels_Araport_geneinfo[,c(1:10,51:58,17:50)]

write.xlsx(KK_GS_arenosa_Araport,"Genes_KletKowa_arenosa_Araportadded.xlsx",sheetName="Classical_genome_scans",col.names=TRUE,row.names=FALSE)
write.xlsx2(KK_GS_arenosa_highSNPs_Araport,"Genes_KletKowa_arenosa_Araportadded.xlsx",sheetName="High_effect_SNPs",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(KK_GS_arenosa_Indels_Araport,"Genes_KletKowa_arenosa_Araportadded.xlsx",sheetName="Indels",col.names=TRUE,row.names=FALSE,append=T)






TH_ortho_numbers<-table(Thal_ortholist$OrthoGroup)
TH_ortho_number<-data.frame(TH_ortho_numbers)
names(TH_ortho_number)<-c("OrthoGroup","Frequency")
Ah_ortho_numbers<-table(Ahal_ortholist$OrthoGroup)
Ah_ortho_number<-data.frame(Ah_ortho_numbers)
names(Ah_ortho_number)<-c("OrthoGroup","Frequency")
Ly_ortho_numbers<-table(Lyr_ortholist$OrthoGroup)
Ly_ortho_number<-data.frame(Ly_ortho_numbers)
names(Ly_ortho_number)<-c("OrthoGroup","Frequency")

data<-read.table("KK_HighSNPs_ortho.csv",sep=";",header=T)
data2<-merge(data,Ah_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)
data3<-merge(data2,TH_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)
data4<-merge(data3,Ly_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)

names(data4)[4:6]<-c("Ah_paralogues","Ath_paralogues","Al_paralogues")
write.table(data4,"KK_higheffectSNPs_paralogues.table",sep="\t",row.names=F)



data<-read.table("MZ_HighSNPs_ortho.csv",sep=";",header=T)
data2<-merge(data,Ah_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)
data3<-merge(data2,TH_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)
data4<-merge(data3,Ly_ortho_number,by.x="OrthoGroup",by.y="OrthoGroup",all.x=T)

names(data4)[4:6]<-c("Ah_paralogues","Ath_paralogues","Al_paralogues")
write.table(data4,"MZ_higheffectSNPs_paralogues.table",sep="\t",row.names=F)



options(java.parameters = "-Xmx24000m")

require(xlsx)

KK_GS_halleri<-read.xlsx("Genes_KletKowa_halleri_Araportadded_W.xlsx",1,header=TRUE)
KK_GS_halleri_highSNPs<-read.xlsx2("Genes_KletKowa_halleri_Araportadded_W.xlsx",2,header=TRUE)
KK_GS_halleri_Indels<-read.xlsx2("Genes_KletKowa_halleri_Araportadded_W.xlsx",3,header=TRUE)

MZ_GS_arenosa<-read.xlsx("Genes_MiasZapa_arenosa_Araportadded_W.xlsx",1,header=TRUE)
MZ_GS_arenosa_highSNPs<-read.xlsx2("Genes_MiasZapa_arenosa_Araportadded_W.xlsx",2,header=TRUE)
MZ_GS_arenosa_Indels<-read.xlsx2("Genes_MiasZapa_arenosa_Araportadded_W.xlsx",3,header=TRUE)

MZ_GS_halleri<-read.xlsx("Genes_MiasZako_halleri_Araportadded_W.xlsx",1,header=TRUE)
MZ_GS_halleri_highSNPs<-read.xlsx2("Genes_MiasZako_halleri_Araportadded_W.xlsx",2,header=TRUE)
MZ_GS_halleri_Indels<-read.xlsx2("Genes_MiasZako_halleri_Araportadded_W.xlsx",3,header=TRUE)

KK_GS_arenosa<-read.xlsx("Genes_KletKowa_arenosa_Araportadded_W.xlsx",1,header=TRUE)
KK_GS_arenosa_highSNPs<-read.xlsx2("Genes_KletKowa_arenosa_Araportadded_W.xlsx",2,header=TRUE)
KK_GS_arenosa_Indels<-read.xlsx2("Genes_KletKowa_arenosa_Araportadded_W.xlsx",3,header=TRUE)
names(KK_GS_halleri_Indels)[32]<-"ID"
names(MZ_GS_arenosa_Indels)[32]<-"ID"
Rank_KK_arenosa<-read.table("Genes_KletKowa_arenosa_rank.txt",sep="\t",header=T)
KK_GS_arenosa_ranked<-merge(KK_GS_arenosa,Rank_KK_arenosa,by.x="AL_ID",by.y="Lyr_Gene",all.x=T)


dup<-duplicated(KK_GS_halleri$AL_ID)
KK_GS_halleri$Rank_in_outlier_metric<-as.numeric(as.character(KK_GS_halleri$Rank_in_outlier_metric))
duplist<-unique(KK_GS_halleri$AL_ID[dup])
KK_GS_halleri_uniquewindows<-KK_GS_halleri[!KK_GS_halleri$AL_ID%in%duplist,]
for (i in 1:length(duplist))
	{Lyrgene<-KK_GS_halleri[KK_GS_halleri$AL_ID==duplist[i],]
	KK_GS_halleri_uniquewindows<-rbind(KK_GS_halleri_uniquewindows,Lyrgene[Lyrgene$Rank_in_outlier_metric==min(Lyrgene$Rank_in_outlier_metric),])
	}
dup<-duplicated(MZ_GS_halleri$AL_ID)
MZ_GS_halleri$Rank_in_outlier_metric<-as.numeric(as.character(MZ_GS_halleri$Rank_in_outlier_metric))
duplist<-unique(MZ_GS_halleri$AL_ID[dup])
MZ_GS_halleri_uniquewindows<-MZ_GS_halleri[!MZ_GS_halleri$AL_ID%in%duplist,]
for (i in 1:length(duplist))
	{Lyrgene<-MZ_GS_halleri[MZ_GS_halleri$AL_ID==duplist[i],]
	MZ_GS_halleri_uniquewindows<-rbind(MZ_GS_halleri_uniquewindows,Lyrgene[Lyrgene$Rank_in_outlier_metric==min(Lyrgene$Rank_in_outlier_metric),])
	}

dup<-duplicated(KK_GS_arenosa_ranked$AL_ID)
KK_GS_arenosa_ranked$Rank_in_outlier_metric<-as.numeric(as.character(KK_GS_arenosa_ranked$Rank_in_outlier_metric))
duplist<-unique(KK_GS_arenosa_ranked$AL_ID[dup])
KK_GS_arenosa_uniquewindows<-KK_GS_arenosa_ranked[!KK_GS_arenosa_ranked$AL_ID%in%duplist,]
for (i in 1:length(duplist))
	{Lyrgene<-KK_GS_arenosa_ranked[KK_GS_arenosa_ranked$AL_ID==duplist[i],]
	KK_GS_arenosa_uniquewindows<-rbind(KK_GS_arenosa_uniquewindows,Lyrgene[Lyrgene$Rank_in_outlier_metric==min(Lyrgene$Rank_in_outlier_metric),])
	}

dup<-duplicated(MZ_GS_arenosa$AL_ID)
MZ_GS_arenosa$Rank_in_outlier_metric<-as.numeric(as.character(MZ_GS_arenosa$Rank_in_outlier_metric))
duplist<-unique(MZ_GS_arenosa$AL_ID[dup])
MZ_GS_arenosa_uniquewindows<-MZ_GS_arenosa[!MZ_GS_arenosa$AL_ID%in%duplist,]
for (i in 1:length(duplist))
	{Lyrgene<-MZ_GS_arenosa[MZ_GS_arenosa$AL_ID==duplist[i],]
	MZ_GS_arenosa_uniquewindows<-rbind(MZ_GS_arenosa_uniquewindows,Lyrgene[Lyrgene$Rank_in_outlier_metric==min(Lyrgene$Rank_in_outlier_metric),])
	}


KK_GS_halleri_highSNPs$Outlier_metric<-rep("High_effect_SNP",nrow(KK_GS_halleri_highSNPs))
KK_GS_halleri_Indels$Outlier_metric<-rep("Indel",nrow(KK_GS_halleri_Indels))
MZ_GS_halleri_highSNPs$Outlier_metric<-rep("High_effect_SNP",nrow(MZ_GS_halleri_highSNPs))
MZ_GS_halleri_Indels$Outlier_metric<-rep("Indel",nrow(MZ_GS_halleri_Indels))
KK_GS_arenosa_highSNPs$Outlier_metric<-rep("High_effect_SNP",nrow(KK_GS_arenosa_highSNPs))
KK_GS_arenosa_Indels$Outlier_metric<-rep("Indel",nrow(KK_GS_arenosa_Indels))
MZ_GS_arenosa_highSNPs$Outlier_metric<-rep("High_effect_SNP",nrow(MZ_GS_arenosa_highSNPs))
MZ_GS_arenosa_Indels$Outlier_metric<-rep("Indel",nrow(MZ_GS_arenosa_Indels))


KK_halleri_MZ_arenosa<-merge(rbind(KK_GS_halleri_uniquewindows[,c(1:10,65:77,63)],KK_GS_halleri_highSNPs[,c(1:10,25:37,41)],KK_GS_halleri_Indels[,c(1:10,32:45)]),rbind(MZ_GS_arenosa_uniquewindows[,c(1:10,61:73,58)],MZ_GS_arenosa_highSNPs[c(1:10,25:37,41)],MZ_GS_arenosa_Indels[,c(1:10,32:45)]),by="AL_ID")
KK_halleri_MZ_arenosa2<-KK_halleri_MZ_arenosa[,c(1:24,47)]
KK_halleri_MZ_arenosa3<-KK_halleri_MZ_arenosa2[!duplicated(KK_halleri_MZ_arenosa2),]

Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)
Ahal_ortholist=read.delim("HalleriGeneOrthogroup.txt", header=TRUE)

TH_ortho_numbers<-table(Thal_ortholist$OrthoGroup)
TH_ortho_number<-data.frame(TH_ortho_numbers)
names(TH_ortho_number)<-c("OrthoGroup","Frequency")
Ah_ortho_numbers<-table(Ahal_ortholist$OrthoGroup)
Ah_ortho_number<-data.frame(Ah_ortho_numbers)
names(Ah_ortho_number)<-c("OrthoGroup","Frequency")
Ly_ortho_numbers<-table(Lyr_ortholist$OrthoGroup)
Ly_ortho_number<-data.frame(Ly_ortho_numbers)
names(Ly_ortho_number)<-c("OrthoGroup","Frequency")

KK_halleri_MZ_arenosa4<-merge(KK_halleri_MZ_arenosa3,Ah_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)
KK_halleri_MZ_arenosa5<-merge(KK_halleri_MZ_arenosa4,TH_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)
KK_halleri_MZ_arenosa6<-merge(KK_halleri_MZ_arenosa5,Ly_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)

names(KK_halleri_MZ_arenosa6)[26:28]<-c("Ah_paralogues","Ath_paralogues","Al_paralogues")


KK_halleri_MZ_arenosa7<-Reduce(function(...) merge(..., by="AL_ID",all.x=TRUE), list(KK_halleri_MZ_arenosa6,KK_GS_halleri_highSNPs[,-c(2:10,25:41)],KK_GS_halleri_Indels[,-c(2:10,32:45)],KK_GS_halleri_uniquewindows[,-c(2:10,63,65:85)],MZ_GS_arenosa_highSNPs[,-c(2:10,25:41)],MZ_GS_arenosa_Indels[,-c(2:10,32:45)],MZ_GS_arenosa_uniquewindows[,-c(2:10,61:81,58)]))
KK_halleri_MZ_arenosa8<-KK_halleri_MZ_arenosa7[!duplicated(KK_halleri_MZ_arenosa7[,c(1:3,24:25,29:200)]),]

write.xlsx(KK_halleri_MZ_arenosa8,"KK_halleri_MZ_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded.xlsx",sheetName="KK_halleri_MZ_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded",row.names=F)











names(KK_GS_arenosa_Indels)[32]<-"ID"
names(MZ_GS_halleri_Indels)[33]<-"ID"






KK_arenosa_MZ_halleri<-merge(rbind(KK_GS_arenosa_uniquewindows[,c(1:10,61:73,59)],KK_GS_arenosa_highSNPs[,c(1:10,25:37,41)],KK_GS_arenosa_Indels[,c(1:10,32:45)]),rbind(MZ_GS_halleri_uniquewindows[,c(1:10,66:78,64)],MZ_GS_halleri_highSNPs[c(1:10,26:38,42)],MZ_GS_halleri_Indels[,c(1:10,33:46)]),by="AL_ID")
KK_arenosa_MZ_halleri2<-KK_arenosa_MZ_halleri[,c(1:24,47)]
KK_arenosa_MZ_halleri3<-KK_arenosa_MZ_halleri2[!duplicated(KK_arenosa_MZ_halleri2),]

Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)
Ahal_ortholist=read.delim("HalleriGeneOrthogroup.txt", header=TRUE)

TH_ortho_numbers<-table(Thal_ortholist$OrthoGroup)
TH_ortho_number<-data.frame(TH_ortho_numbers)
names(TH_ortho_number)<-c("OrthoGroup","Frequency")
Ah_ortho_numbers<-table(Ahal_ortholist$OrthoGroup)
Ah_ortho_number<-data.frame(Ah_ortho_numbers)
names(Ah_ortho_number)<-c("OrthoGroup","Frequency")
Ly_ortho_numbers<-table(Lyr_ortholist$OrthoGroup)
Ly_ortho_number<-data.frame(Ly_ortho_numbers)
names(Ly_ortho_number)<-c("OrthoGroup","Frequency")

KK_arenosa_MZ_halleri4<-merge(KK_arenosa_MZ_halleri3,Ah_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)
KK_arenosa_MZ_halleri5<-merge(KK_arenosa_MZ_halleri4,TH_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)
KK_arenosa_MZ_halleri6<-merge(KK_arenosa_MZ_halleri5,Ly_ortho_number,by.x="OrthoGroup.x",by.y="OrthoGroup",all.x=T)

names(KK_arenosa_MZ_halleri6)[26:28]<-c("Ah_paralogues","Ath_paralogues","Al_paralogues")


KK_arenosa_MZ_halleri7<-Reduce(function(...) merge(..., by="AL_ID",all.x=TRUE), list(KK_arenosa_MZ_halleri6,KK_GS_arenosa_highSNPs[,-c(2:10,25:41)],KK_GS_arenosa_Indels[,-c(2:10,32:45)],KK_GS_arenosa_uniquewindows[,-c(2:10,59,61:82)],MZ_GS_halleri_highSNPs[,-c(2:10,26:42)],MZ_GS_halleri_Indels[,-c(2:10,33:46)],MZ_GS_halleri_uniquewindows[,-c(2:10,66:86,64)]))
KK_arenosa_MZ_halleri8<-KK_arenosa_MZ_halleri7[!duplicated(KK_arenosa_MZ_halleri7[,c(1:3,24:25,29:203)]),]



geneinfo<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",header=F)
names(geneinfo)<-c("Scaffold","Software","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
geneinfo<-geneinfo[geneinfo$Type=="gene",]
ID<-strsplit(as.character(geneinfo$ID),"=",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T))
geneinfo$ID<-as.character(ID2$X3)
geneinfo$ID<-as.character(geneinfo$ID)
ID<-strsplit(as.character(geneinfo$ID),";",fixed=T)
ID1<-lapply(ID,"[",1:max(unlist(lapply(ID,length))))
ID2<-data.frame(matrix(unlist(ID1),nrow=length(ID1),byrow=T))
geneinfo$ID<-as.character(ID2$X1)

KK_arenosa_MZ_halleri8_geneinfo<-merge(KK_arenosa_MZ_halleri8,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_arenosa_MZ_halleri8_geneinfo2<-KK_arenosa_MZ_halleri8_geneinfo[,c(1:10,204:211,11:203)]
write.xlsx(KK_arenosa_MZ_halleri8_geneinfo2,"KK_arenosa_MZ_halleri_overlap_lyratagenes_filtered0_1perc_Araportadded.xlsx",sheetName="KK_arenosa_MZ_halleri_overlap_lyratagenes_filtered0_1perc_Araportadded",row.names=F)

KK_halleri_MZ_arenosa8_geneinfo<-merge(KK_halleri_MZ_arenosa8,geneinfo,by.x="AL_ID",by.y="ID",no.dups=F,all.x=T)
KK_halleri_MZ_arenosa8_geneinfo2<-KK_halleri_MZ_arenosa8_geneinfo[,c(1:10,201:208,11:200)]
write.xlsx(KK_halleri_MZ_arenosa8_geneinfo2,"KK_halleri_MZ_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded.xlsx",sheetName="KK_halleri_MZ_arenosa_overlap_lyratagenes_filtered0_1perc_Araportadded",row.names=F)

