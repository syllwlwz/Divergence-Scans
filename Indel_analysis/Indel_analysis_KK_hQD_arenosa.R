data<-read.csv("LyrataKletKowaGS_INDEL_hQD_arenosa.csv",sep="\t",header=T)
KK_Indels<-data[!((data$AC/36)==(data$AC.1/32)),]

SNPeff<-read.csv("KletKowa_INDEL_hQD_ann.arenosa.table",header=T,sep="\t")
KK_HC<-read.table("KletKowa_INDEL_hQD.arenosa.ann.vcf",sep="\t",header=F)

SNPeff2<-data.frame(SNPeff,KK_HC[,4:5])
names(SNPeff2)[4:5]<-c("REF","ALT")

SNPeffAnnsplit<-strsplit(as.character(SNPeff2$ANN),",",fixed=T)
SNPeffAnnsplit1<-lapply(SNPeffAnnsplit,"[",1:max(unlist(lapply(SNPeffAnnsplit,length))))
SNPeffAnnsplit2<-data.frame(matrix(unlist(SNPeffAnnsplit1),nrow=length(SNPeffAnnsplit1),byrow=T))
SNPeffdf<-data.frame()
for (i in 1:max(unlist(lapply(SNPeffAnnsplit,length))))
	{SNPeff1<-ifelse(!is.na(SNPeffAnnsplit2[,i]),strsplit(as.character(SNPeffAnnsplit2[,i]),"|",fixed=T),NA)
	SNPeffAnn1<-lapply(SNPeff1,"[",1:4)
	SNPeffAnndf1<-data.frame(matrix(unlist(SNPeffAnn1),nrow=length(SNPeffAnnsplit1),byrow=T))
	names(SNPeffAnndf1)<-c("Base","Characterization","Effect","Gene")
	SNPeffcomp1<-cbind(SNPeff2[,c(1:2,4:5)],SNPeffAnndf1)
	SNPeffdf<-rbind(SNPeffdf,na.exclude(SNPeffcomp1))
	}
require(GenomicRanges)

SNPeffdf_GRange<-GRanges(seqnames=SNPeffdf$CHROM,ranges=IRanges(start=SNPeffdf$POS,end=SNPeffdf$POS))
values(SNPeffdf_GRange)<-SNPeffdf[,3:8]

KK_Indels_GRange<-GRanges(seqnames=KK_Indels$Chrom,ranges=IRanges(start=KK_Indels$POS,end=KK_Indels$POS))
values(KK_Indels_GRange)<-KK_Indels[,c("AC","AC.1")]
SNP_KK_Indels=mergeByOverlaps(KK_Indels_GRange,SNPeffdf_GRange,type=c("any"))
SNP_KK_Indels_df=as.data.frame(SNP_KK_Indels)
SNP_KK_Indels_df2<-cbind(SNP_KK_Indels_df[,1:2],SNP_KK_Indels_df[,8:9],SNP_KK_Indels_df[,21:26])
names(SNP_KK_Indels_df2)<-c("Chrom","Pos","AC_Klet","AC_Kowa","REF","ALT","Base","Characterization","Effect","Gene")

SNP_KK_Indels_diff<-SNP_KK_Indels_df2[abs((SNP_KK_Indels_df2$AC_Klet/36)-(SNP_KK_Indels_df2$AC_Kowa/32))>0.4&(((SNP_KK_Indels_df2$AC_Klet/36)<=0.1)|((SNP_KK_Indels_df2$AC_Klet/36)>=0.9)),]

Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)
#Ahal_ortholist=read.delim("HalleriGeneOrthogroup.txt", header=TRUE)
Thal_description=read.delim("TAIR10_functional_descriptions.txt", header=TRUE)
colnames(Thal_description)[1]<-"Gene"
test1<-substr(Thal_description[,1],1,9)
Thal_description[,1]<-test1

SNP_KK_Indels_diff_w_orthogroup= merge(SNP_KK_Indels_diff, Lyr_ortholist, by="Gene", all.x=TRUE)
SNP_KK_Indels_diff_w_orthogroup_Thalgenes= merge(SNP_KK_Indels_diff_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(SNP_KK_Indels_diff_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(SNP_KK_Indels_diff_w_orthogroup_Thalgenes)[12]="Gene"
Thalgenes_described_SNP_KK_Indels_diff_w_orthogroup= merge(Thal_description, SNP_KK_Indels_diff_w_orthogroup_Thalgenes, by ="Gene",all.y=T)

require(xlsx)
Methomlist<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",2,header=T)

Thalgenes_described_SNP_KK_Indels_diff_w_orthogroup$Gene<-toupper(Thalgenes_described_SNP_KK_Indels_diff_w_orthogroup$Gene)
Methomlist$AGI_code<-toupper(Methomlist$AGI_code)
Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist<-merge(Thalgenes_described_SNP_KK_Indels_diff_w_orthogroup,Methomlist,all.x=T,by.x="Gene",by.y="AGI_code")
Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist$AC_Klet<-Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist$AC_Klet/36
Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist$AC_Kowa<-Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist$AC_Kowa/32
names(Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist)[10]<-"AF_Klet"
names(Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist)[11]<-"AF_Kowa"
write.table(Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist,"SNP_KK_Indels_diffabove04_hQD_arenosa.txt", sep="\t", row.names=F)

#GROM
Klet_indel_grom<-read.table("KletKowa_INDEL_GROM.ann.vcf",sep="\t",header=F)
Kowa_indel_grom<-read.table("KletKowa_INDEL_GROM.Kowa.ann.vcf",sep="\t",header=F)
names(Klet_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","ANN","FORMAT_names","FORMAT")
names(Kowa_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","ANN","FORMAT_names","FORMAT")

Klet_indel_gromAnnsplit<-strsplit(as.character(Klet_indel_grom$ANN),",",fixed=T)
Klet_indel_gromAnnsplit1<-lapply(Klet_indel_gromAnnsplit,"[",1:max(unlist(lapply(Klet_indel_gromAnnsplit,length))))
Klet_indel_gromAnnsplit2<-data.frame(matrix(unlist(Klet_indel_gromAnnsplit1),nrow=length(Klet_indel_gromAnnsplit1),byrow=T))
Klet_indel_gromdf<-data.frame()
for (i in 1:max(unlist(lapply(Klet_indel_gromAnnsplit,length))))
	{Klet_indel_grom1<-ifelse(!is.na(Klet_indel_gromAnnsplit2[,i]),strsplit(as.character(Klet_indel_gromAnnsplit2[,i]),"|",fixed=T),NA)
	Klet_indel_gromAnn1<-lapply(Klet_indel_grom1,"[",1:4)
	Klet_indel_gromAnndf1<-data.frame(matrix(unlist(Klet_indel_gromAnn1),nrow=length(Klet_indel_gromAnnsplit1),byrow=T))
	names(Klet_indel_gromAnndf1)<-c("Base","Characterization","Effect","Gene")
	Klet_indel_gromcomp1<-cbind(Klet_indel_grom[,1:2],Klet_indel_gromAnndf1)
	Klet_indel_gromdf<-rbind(Klet_indel_gromdf,na.exclude(Klet_indel_gromcomp1))
	}
require(GenomicRanges)

Klet_indel_gromdfMODIFIER<-Klet_indel_gromdf[Klet_indel_gromdf$Effect=="MODIFIER",]
Klet_indel_gromdfMODIFIER_GRange<-GRanges(seqnames=Klet_indel_gromdfMODIFIER$CHROM,ranges=IRanges(start=Klet_indel_gromdfMODIFIER$POS,end=Klet_indel_gromdfMODIFIER$POS))
values(Klet_indel_gromdfMODIFIER_GRange)<-Klet_indel_gromdfMODIFIER[,1:6]

Klet_indel_gromdfMODERATE<-Klet_indel_gromdf[Klet_indel_gromdf$Effect=="MODERATE",]
Klet_indel_gromdfMODERATE_GRange<-GRanges(seqnames=Klet_indel_gromdfMODERATE$CHROM,ranges=IRanges(start=Klet_indel_gromdfMODERATE$POS,end=Klet_indel_gromdfMODERATE$POS))
values(Klet_indel_gromdfMODERATE_GRange)<-Klet_indel_gromdfMODERATE[,1:6]

Klet_indel_gromdfHIGH<-Klet_indel_gromdf[Klet_indel_gromdf$Effect=="HIGH",]
Klet_indel_gromdfHIGH_GRange<-GRanges(seqnames=Klet_indel_gromdfHIGH$CHROM,ranges=IRanges(start=Klet_indel_gromdfHIGH$POS,end=Klet_indel_gromdfHIGH$POS))
values(Klet_indel_gromdfHIGH_GRange)<-Klet_indel_gromdfHIGH[,3:6]

KK_Indels_GRange<-GRanges(seqnames=KK_Indels$Chrom,ranges=IRanges(start=KK_Indels$POS,end=KK_Indels$POS))
values(KK_Indels_GRange)<-KK_Indels[,c("AC","AC.1")]
SNPHIGH_KK_Indels=mergeByOverlaps(KK_Indels_GRange,SNPeffdfHIGH_GRange,type=c("any"))
SNPHIGH_KK_Indels_df=as.data.frame(SNPHIGH_KK_Indels)
SNPHIGH_KK_Indels_df2<-cbind(SNPHIGH_KK_Indels_df[,1:2],SNPHIGH_KK_Indels_df[,8:9],SNPHIGH_KK_Indels_df[,19:22])
names(SNPHIGH_KK_Indels_df2)<-c("Chrom","Pos","AC_Klet","AC_Kowa","Base","Characterization","Effect","Gene")

SNPHIGH_KK_Indels_diff<-SNPHIGH_KK_Indels_df2[abs((SNPHIGH_KK_Indels_df2$AC_Klet/20)-(SNPHIGH_KK_Indels_df2$AC_Kowa/14))>0.4,]

Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)
Ahal_ortholist=read.delim("HalleriGeneOrthogroup.txt", header=TRUE)
Thal_description=read.delim("TAIR10_functional_descriptions.txt", header=TRUE)
colnames(Thal_description)[1]<-"Gene"
test1<-substr(Thal_description[,1],1,9)
Thal_description[,1]<-test1

SNPHIGH_KK_Indels_diff_w_orthogroup= merge(SNPHIGH_KK_Indels_diff, Lyr_ortholist, by="Gene", all.x=TRUE)
SNPHIGH_KK_Indels_diff_w_orthogroup_Thalgenes= merge(SNPHIGH_KK_Indels_diff_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(SNPHIGH_KK_Indels_diff_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(SNPHIGH_KK_Indels_diff_w_orthogroup_Thalgenes)[10]="Gene"
Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroup= merge(Thal_description, SNPHIGH_KK_Indels_diff_w_orthogroup_Thalgenes, by ="Gene",all.y=T)

require(xlsx)
Methomlist<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",2,header=T)

Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroup$Gene<-toupper(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroup$Gene)
Methomlist$AGI_code<-toupper(Methomlist$AGI_code)
Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist<-merge(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroup,Methomlist,all.x=T,by.x="Gene",by.y="AGI_code")
Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$AC_Klet<-Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$AC_Klet/20
Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$AC_Kowa<-Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$AC_Kowa/14
names(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist)[10]<-"AF_Klet"
names(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist)[11]<-"AF_Kowa"
write.table(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist,"SNPHIGH_KK_Indels_diffabove04_hQD.txt", sep="\t", row.names=F)

Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist_fix<-Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist[Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$AF_Klet>=0.9|Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$AF_Klet<=0.1,]

write.table(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist_fix,"SNPHIGH_KK_Indels_diffabove04_fixed_hQD.txt", sep="\t", row.names=F)

length(unique(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist_fix$Gene))
#1272

SNPeffdf_GRange<-GRanges(seqnames=SNPeffdf$CHROM,ranges=IRanges(start=SNPeffdf$POS,end=SNPeffdf$POS))
values(SNPeffdf_GRange)<-SNPeffdf[,1:6]

SNPHIGH_KK_Indels=mergeByOverlaps(KK_Indels_GRange,SNPeffdf_GRange,type=c("any"))
SNPHIGH_KK_Indels_df=as.data.frame(SNPHIGH_KK_Indels)
SNPHIGH_KK_Indels_df2<-cbind(SNPHIGH_KK_Indels_df[,1:2],SNPHIGH_KK_Indels_df[,8:9],SNPHIGH_KK_Indels_df[,23:26])
names(SNPHIGH_KK_Indels_df2)<-c("Chrom","Pos","AC_Klet","AC_Kowa","Base","Characterization","Effect","Gene")

SNPHIGH_KK_Indels_diff<-SNPHIGH_KK_Indels_df2[abs((SNPHIGH_KK_Indels_df2$AC_Klet/20)-(SNPHIGH_KK_Indels_df2$AC_Kowa/14))>0.4,]

Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)
Ahal_ortholist=read.delim("HalleriGeneOrthogroup.txt", header=TRUE)
Thal_description=read.delim("TAIR10_functional_descriptions.txt", header=TRUE)
colnames(Thal_description)[1]<-"Gene"
test1<-substr(Thal_description[,1],1,9)
Thal_description[,1]<-test1

SNPHIGH_KK_Indels_diff_w_orthogroup= merge(SNPHIGH_KK_Indels_diff, Lyr_ortholist, by="Gene", all.x=TRUE)
SNPHIGH_KK_Indels_diff_w_orthogroup_Thalgenes= merge(SNPHIGH_KK_Indels_diff_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(SNPHIGH_KK_Indels_diff_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(SNPHIGH_KK_Indels_diff_w_orthogroup_Thalgenes)[10]="Gene"
Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroup= merge(Thal_description, SNPHIGH_KK_Indels_diff_w_orthogroup_Thalgenes, by ="Gene",all.y=T)

require(xlsx)
Methomlist<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",2,header=T)

Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroup$Gene<-toupper(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroup$Gene)
Methomlist$AGI_code<-toupper(Methomlist$AGI_code)
Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist<-merge(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroup,Methomlist,all.x=T,by.x="Gene",by.y="AGI_code")
Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$AC_Klet<-Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$AC_Klet/20
Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$AC_Kowa<-Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$AC_Kowa/14
names(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist)[10]<-"AF_Klet"
names(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist)[11]<-"AF_Kowa"
Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist2<-Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist[!Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist$Effect=="LOW",]
write.table(Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist,"SNPHIGH_KK_Indels_diffabove04_hQD.txt", sep="\t", row.names=F)

test<-Thalgenes_described_SNPHIGH_KK_Indels_diff_w_orthogroupmetlist2
test_GR<-GRanges(seqnames=tolower(test$Chrom),ranges=IRanges(start=test$Pos,end=test$Pos))
test_nearest<-nearest(test_GR,KletKowa_Indels_GR,ignore.strand=T)
test3<-KletKowa_Indels[test_nearest,1:10]
test4<-cbind(test,test3)
test5<-test4[abs(test4[,9]-test4[,29])<=nchar(as.character(test4$Base)),]

