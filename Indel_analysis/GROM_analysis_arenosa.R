##########################################################
#Load and format coverage data per gene for lyrata genome#
##########################################################
options(java.parameters = "-Xmx20000m")

Mias_genes_cov<-read.table("Hist_arenosa_Mias.bed",sep="\t",quote="",fill=T)
Mias_genes_cov<-Mias_genes_cov[Mias_genes_cov$V1!="all",]
names(Mias_genes_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Gene_width","Fraction_of_gene_with_coverage")
Mias_genes_cov2<-Mias_genes_cov[Mias_genes_cov$Type=="gene",]
Mias_genes_cov2<-droplevels(Mias_genes_cov2)
Mias_genes_cov2$ID<-gsub("ID=","",gsub(";.*","",Mias_genes_cov2$ID))
mean_cov<-function(arg1,arg2){
	meancov<-mean(rep(arg1,arg2))
	return(meancov)
}
median_cov<-function(arg1,arg2){
	mediancov<-as.double(median(rep(arg1,arg2)))
	return(mediancov)
}

require(data.table)
Mias_genes_cov_dt<-data.table(Mias_genes_cov2)
Mias_genes_cov_dt<-droplevels(Mias_genes_cov_dt)
Mias_mean_cov<-Mias_genes_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by="ID"]
Mias_median_cov<-Mias_genes_cov_dt[,median_cov(Coverage,Bases_at_coverage),by="ID"]
Mias_genes_coverage<-data.frame(unique(Mias_genes_cov2[,1:9]),Mias_mean_cov[,2],Mias_median_cov[,2])
names(Mias_genes_coverage)[10:11]<-c("Mean_coverage_gene","Median_coverage_gene")

Zapa_genes_cov<-read.table("Hist_arenosa_Zapa.bed",sep="\t",quote="",fill=T)
Zapa_genes_cov<-Zapa_genes_cov[Zapa_genes_cov$V1!="all",]
names(Zapa_genes_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Gene_width","Fraction_of_gene_with_coverage")
Zapa_genes_cov2<-Zapa_genes_cov[Zapa_genes_cov$Type=="gene",]
Zapa_genes_cov2<-droplevels(Zapa_genes_cov2)
Zapa_genes_cov2$ID<-gsub("ID=","",gsub(";.*","",Zapa_genes_cov2$ID))
Zapa_genes_cov_dt<-data.table(Zapa_genes_cov2)
Zapa_genes_cov_dt<-droplevels(Zapa_genes_cov_dt)
Zapa_mean_cov<-Zapa_genes_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by="ID"]
Zapa_median_cov<-Zapa_genes_cov_dt[,median_cov(Coverage,Bases_at_coverage),by="ID"]
Zapa_genes_coverage<-data.frame(unique(Zapa_genes_cov2[,1:9]),Zapa_mean_cov[,2],Zapa_median_cov[,2])
names(Zapa_genes_coverage)[10:11]<-c("Mean_coverage_gene","Median_coverage_gene")

Klet_genes_cov<-read.table("Hist_arenosa_Klet.bed",sep="\t",quote="",fill=T)
Klet_genes_cov<-Klet_genes_cov[Klet_genes_cov$V1!="all",]
names(Klet_genes_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Gene_width","Fraction_of_gene_with_coverage")
Klet_genes_cov2<-Klet_genes_cov[Klet_genes_cov$Type=="gene",]
Klet_genes_cov2<-droplevels(Klet_genes_cov2)
Klet_genes_cov2$ID<-gsub("ID=","",gsub(";.*","",Klet_genes_cov2$ID))
Klet_genes_cov_dt<-data.table(Klet_genes_cov2)
Klet_genes_cov_dt<-droplevels(Klet_genes_cov_dt)
Klet_mean_cov<-Klet_genes_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by="ID"]
Klet_median_cov<-Klet_genes_cov_dt[,median_cov(Coverage,Bases_at_coverage),by="ID"]
Klet_genes_coverage<-data.frame(unique(Klet_genes_cov2[,1:9]),Klet_mean_cov[,2],Klet_median_cov[,2])
names(Klet_genes_coverage)[10:11]<-c("Mean_coverage_gene","Median_coverage_gene")

Kowa_genes_cov<-read.table("Hist_arenosa_Kowa.bed",sep="\t",quote="",fill=T)
Kowa_genes_cov<-Kowa_genes_cov[Kowa_genes_cov$V1!="all",]
names(Kowa_genes_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Gene_width","Fraction_of_gene_with_coverage")
Kowa_genes_cov2<-Kowa_genes_cov[Kowa_genes_cov$Type=="gene",]
Kowa_genes_cov2<-droplevels(Kowa_genes_cov2)
Kowa_genes_cov2$ID<-gsub("ID=","",gsub(";.*","",Kowa_genes_cov2$ID))
Kowa_genes_cov_dt<-data.table(Kowa_genes_cov2)
Kowa_genes_cov_dt<-droplevels(Kowa_genes_cov_dt)
Kowa_mean_cov<-Kowa_genes_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by="ID"]
Kowa_median_cov<-Kowa_genes_cov_dt[,median_cov(Coverage,Bases_at_coverage),by="ID"]
Kowa_genes_coverage<-data.frame(unique(Kowa_genes_cov2[,1:9]),Kowa_mean_cov[,2],Kowa_median_cov[,2])
names(Kowa_genes_coverage)[10:11]<-c("Mean_coverage_gene","Median_coverage_gene")

Wulm_genes_cov<-read.table("Hist_arenosa_Wulm.bed",sep="\t",quote="",fill=T)
Wulm_genes_cov<-Wulm_genes_cov[Wulm_genes_cov$V1!="all",]
names(Wulm_genes_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Gene_width","Fraction_of_gene_with_coverage")
Wulm_genes_cov2<-Wulm_genes_cov[Wulm_genes_cov$Type=="gene",]
Wulm_genes_cov2<-droplevels(Wulm_genes_cov2)
Wulm_genes_cov2$ID<-gsub("ID=","",gsub(";.*","",Wulm_genes_cov2$ID))
Wulm_genes_cov_dt<-data.table(Wulm_genes_cov2)
Wulm_genes_cov_dt<-droplevels(Wulm_genes_cov_dt)
Wulm_mean_cov<-Wulm_genes_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by="ID"]
Wulm_median_cov<-Wulm_genes_cov_dt[,median_cov(Coverage,Bases_at_coverage),by="ID"]
Wulm_genes_coverage<-data.frame(unique(Wulm_genes_cov2[,1:9]),Wulm_mean_cov[,2],Wulm_median_cov[,2])
names(Wulm_genes_coverage)[10:11]<-c("Mean_coverage_gene","Median_coverage_gene")

Chok_genes_cov<-read.table("Hist_arenosa_CHO.bed",sep="\t",quote="",fill=T)
Chok_genes_cov<-Chok_genes_cov[Chok_genes_cov$V1!="all",]
names(Chok_genes_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Gene_width","Fraction_of_gene_with_coverage")
Chok_genes_cov2<-Chok_genes_cov[Chok_genes_cov$Type=="gene",]
Chok_genes_cov2<-droplevels(Chok_genes_cov2)
Chok_genes_cov2$ID<-gsub("ID=","",gsub(";.*","",Chok_genes_cov2$ID))
Chok_genes_cov_dt<-data.table(Chok_genes_cov2)
Chok_genes_cov_dt<-droplevels(Chok_genes_cov_dt)
Chok_mean_cov<-Chok_genes_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by="ID"]
Chok_median_cov<-Chok_genes_cov_dt[,median_cov(Coverage,Bases_at_coverage),by="ID"]
Chok_genes_coverage<-data.frame(unique(Chok_genes_cov2[,1:9]),Chok_mean_cov[,2],Chok_median_cov[,2])
names(Chok_genes_coverage)[10:11]<-c("Mean_coverage_gene","Median_coverage_gene")


####################################################################
#Load and format coverage data per exons of genes for lyrata genome#
####################################################################
Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)
Ahal_ortholist=read.delim("HalleriGeneOrthogroup.txt", header=TRUE)
Thal_description=read.delim("TAIR10_functional_descriptions.txt", header=TRUE)
colnames(Thal_description)[1]<-"Gene"
test1<-substr(Thal_description[,1],1,9)
Thal_description[,1]<-test1
require(xlsx)
Methomlist<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",2,header=T)
Methomlist[,2]<-toupper(as.character(Methomlist[,2]))

Mias_exons_cov<-read.table("Hist_arenosa_exons_Mias.bed",sep="\t",quote="",fill=T)
Mias_exons_cov<-Mias_exons_cov[Mias_exons_cov$V1!="all",]
names(Mias_exons_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Exon_width","Fraction_of_exon_with_coverage")
Mias_exons_cov2<-Mias_exons_cov[Mias_exons_cov$Type=="exon",]
Mias_exons_cov2<-droplevels(Mias_exons_cov2)
Mias_exons_cov2$ID<-gsub("\\.t.","",gsub("ID=","",gsub(";.*","",Mias_exons_cov2$ID)))

mean_cov<-function(arg1,arg2){
	meancov<-mean(rep(arg1,arg2))
	return(meancov)
}
median_cov<-function(arg1,arg2){
	mediancov<-as.double(median(rep(arg1,arg2)))
	return(mediancov)
}

require(data.table)
Mias_exons_cov_dt<-data.table(Mias_exons_cov2)
Mias_exons_cov_dt<-droplevels(Mias_exons_cov_dt)
Mias_mean_exon_cov<-Mias_exons_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Mias_median_exon_cov<-Mias_exons_cov_dt[,median_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Mias_mean_exon_mean_cov<-Mias_mean_exon_cov[,mean_cov(V1,Exon_width),by="ID"]
Mias_median_exon_mean_cov<-Mias_median_exon_cov[,median_cov(V1,Exon_width),by="ID"]

Mias_exons_coverage<-merge(merge(Mias_genes_coverage,Mias_mean_exon_mean_cov,by="ID"),Mias_median_exon_mean_cov,by="ID")
names(Mias_exons_coverage)[12:13]<-c("Mean_coverage_exons","Median_coverage_exons")
Mias_exons_coverage_w_orthogroup= merge(Mias_exons_coverage, Lyr_ortholist, by.x="ID",by.y="Gene", all.x=TRUE)
Mias_exons_coverage_w_orthogroup_Thalgenes= merge(Mias_exons_coverage_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Mias_exons_coverage_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Mias_exons_coverage_w_orthogroup_Thalgenes)[15]="Gene"
Thalgenes_described_Mias_exons_coverage_w_orthogroup= merge(Thal_description, Mias_exons_coverage_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Mias_coverage_all_genes<-merge(Thalgenes_described_Mias_exons_coverage_w_orthogroup,Methomlist,by.x="Gene",by.y="AGI_code",all.x=TRUE)
write.table(Mias_coverage_all_genes,"Mias_coverage_genes_exons_metalllist.table",sep="\t",row.names=F)
write.table(Mias_exons_coverage,"Mias_arenosa_genes_exons.table",sep="\t",row.names=F)

Zapa_exons_cov<-read.table("Hist_arenosa_exons_Zapa.bed",sep="\t",quote="",fill=T)
Zapa_exons_cov<-Zapa_exons_cov[Zapa_exons_cov$V1!="all",]
names(Zapa_exons_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Exon_width","Fraction_of_exon_with_coverage")
Zapa_exons_cov2<-Zapa_exons_cov[Zapa_exons_cov$Type=="exon",]
Zapa_exons_cov2<-droplevels(Zapa_exons_cov2)
Zapa_exons_cov2$ID<-gsub("\\.t.","",gsub("ID=","",gsub(";.*","",Zapa_exons_cov2$ID)))
Zapa_exons_cov_dt<-data.table(Zapa_exons_cov2)
Zapa_exons_cov_dt<-droplevels(Zapa_exons_cov_dt)
Zapa_mean_exon_cov<-Zapa_exons_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Zapa_median_exon_cov<-Zapa_exons_cov_dt[,median_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Zapa_mean_exon_mean_cov<-Zapa_mean_exon_cov[,mean_cov(V1,Exon_width),by="ID"]
Zapa_median_exon_mean_cov<-Zapa_median_exon_cov[,median_cov(V1,Exon_width),by="ID"]
Zapa_exons_coverage<-merge(merge(Zapa_genes_coverage,Zapa_mean_exon_mean_cov,by="ID"),Zapa_median_exon_mean_cov,by="ID")
names(Zapa_exons_coverage)[12:13]<-c("Mean_coverage_exons","Median_coverage_exons")
write.table(Zapa_exons_coverage,"Zapa_arenosa_genes_exons.table",sep="\t",row.names=F)

Klet_exons_cov<-read.table("Hist_arenosa_exons_Klet.bed",sep="\t",quote="",fill=T)
Klet_exons_cov<-Klet_exons_cov[Klet_exons_cov$V1!="all",]
names(Klet_exons_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Exon_width","Fraction_of_exon_with_coverage")
Klet_exons_cov2<-Klet_exons_cov[Klet_exons_cov$Type=="exon",]
Klet_exons_cov2<-droplevels(Klet_exons_cov2)
Klet_exons_cov2$ID<-gsub("\\.t.","",gsub("ID=","",gsub(";.*","",Klet_exons_cov2$ID)))
Klet_exons_cov_dt<-data.table(Klet_exons_cov2)
Klet_exons_cov_dt<-droplevels(Klet_exons_cov_dt)
Klet_mean_exon_cov<-Klet_exons_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Klet_median_exon_cov<-Klet_exons_cov_dt[,median_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Klet_mean_exon_mean_cov<-Klet_mean_exon_cov[,mean_cov(V1,Exon_width),by="ID"]
Klet_median_exon_mean_cov<-Klet_median_exon_cov[,median_cov(V1,Exon_width),by="ID"]
Klet_exons_coverage<-merge(merge(Klet_genes_coverage,Klet_mean_exon_mean_cov,by="ID"),Klet_median_exon_mean_cov,by="ID")
names(Klet_exons_coverage)[12:13]<-c("Mean_coverage_exons","Median_coverage_exons")
write.table(Klet_exons_coverage,"Klet_arenosa_genes_exons.table",sep="\t",row.names=F)

Kowa_exons_cov<-read.table("Hist_arenosa_exons_Kowa.bed",sep="\t",quote="",fill=T)
Kowa_exons_cov<-Kowa_exons_cov[Kowa_exons_cov$V1!="all",]
names(Kowa_exons_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Exon_width","Fraction_of_exon_with_coverage")
Kowa_exons_cov2<-Kowa_exons_cov[Kowa_exons_cov$Type=="exon",]
Kowa_exons_cov2<-droplevels(Kowa_exons_cov2)
Kowa_exons_cov2$ID<-gsub("\\.t.","",gsub("ID=","",gsub(";.*","",Kowa_exons_cov2$ID)))
Kowa_exons_cov_dt<-data.table(Kowa_exons_cov2)
Kowa_exons_cov_dt<-droplevels(Kowa_exons_cov_dt)
Kowa_mean_exon_cov<-Kowa_exons_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Kowa_median_exon_cov<-Kowa_exons_cov_dt[,median_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Kowa_mean_exon_mean_cov<-Kowa_mean_exon_cov[,mean_cov(V1,Exon_width),by="ID"]
Kowa_median_exon_mean_cov<-Kowa_median_exon_cov[,median_cov(V1,Exon_width),by="ID"]
Kowa_exons_coverage<-merge(merge(Kowa_genes_coverage,Kowa_mean_exon_mean_cov,by="ID"),Kowa_median_exon_mean_cov,by="ID")
names(Kowa_exons_coverage)[12:13]<-c("Mean_coverage_exons","Median_coverage_exons")
write.table(Kowa_exons_coverage,"Kowa_arenosa_genes_exons.table",sep="\t",row.names=F)

Wulm_exons_cov<-read.table("Hist_arenosa_exons_Wulm.bed",sep="\t",quote="",fill=T)
Wulm_exons_cov<-Wulm_exons_cov[Wulm_exons_cov$V1!="all",]
names(Wulm_exons_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Exon_width","Fraction_of_exon_with_coverage")
Wulm_exons_cov2<-Wulm_exons_cov[Wulm_exons_cov$Type=="exon",]
Wulm_exons_cov2<-droplevels(Wulm_exons_cov2)
Wulm_exons_cov2$ID<-gsub("\\.t.","",gsub("ID=","",gsub(";.*","",Wulm_exons_cov2$ID)))
Wulm_exons_cov_dt<-data.table(Wulm_exons_cov2)
Wulm_exons_cov_dt<-droplevels(Wulm_exons_cov_dt)
Wulm_mean_exon_cov<-Wulm_exons_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Wulm_median_exon_cov<-Wulm_exons_cov_dt[,median_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Wulm_mean_exon_mean_cov<-Wulm_mean_exon_cov[,mean_cov(V1,Exon_width),by="ID"]
Wulm_median_exon_mean_cov<-Wulm_median_exon_cov[,median_cov(V1,Exon_width),by="ID"]
Wulm_exons_coverage<-merge(merge(Wulm_genes_coverage,Wulm_mean_exon_mean_cov,by="ID"),Wulm_median_exon_mean_cov,by="ID")
names(Wulm_exons_coverage)[12:13]<-c("Mean_coverage_exons","Median_coverage_exons")
write.table(Wulm_exons_coverage,"Wulm_arenosa_genes_exons.table",sep="\t",row.names=F)

Chok_exons_cov<-read.table("Hist_arenosa_exons_CHO.bed",sep="\t",quote="",fill=T)
Chok_exons_cov<-Chok_exons_cov[Chok_exons_cov$V1!="all",]
names(Chok_exons_cov)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID","Coverage","Bases_at_coverage","Exon_width","Fraction_of_exon_with_coverage")
Chok_exons_cov2<-Chok_exons_cov[Chok_exons_cov$Type=="exon",]
Chok_exons_cov2<-droplevels(Chok_exons_cov2)
Chok_exons_cov2$ID<-gsub("\\.t.","",gsub("ID=","",gsub(";.*","",Chok_exons_cov2$ID)))
Chok_exons_cov_dt<-data.table(Chok_exons_cov2)
Chok_exons_cov_dt<-droplevels(Chok_exons_cov_dt)
Chok_mean_exon_cov<-Chok_exons_cov_dt[,mean_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Chok_median_exon_cov<-Chok_exons_cov_dt[,median_cov(Coverage,Bases_at_coverage),by=list(ID,Start_pos,Exon_width)]
Chok_mean_exon_mean_cov<-Chok_mean_exon_cov[,mean_cov(V1,Exon_width),by="ID"]
Chok_median_exon_mean_cov<-Chok_median_exon_cov[,median_cov(V1,Exon_width),by="ID"]
Chok_exons_coverage<-merge(merge(Chok_genes_coverage,Chok_mean_exon_mean_cov,by="ID"),Chok_median_exon_mean_cov,by="ID")
names(Chok_exons_coverage)[12:13]<-c("Mean_coverage_exons","Median_coverage_exons")
write.table(Chok_exons_coverage,"Chok_arenosa_genes_exons.table",sep="\t",row.names=F)


#########################################################################################
Mias_coverage<-read.table("Mias_arenosa_genes_exons.table",sep="\t",header=T)
Zapa_coverage<-read.table("Zapa_arenosa_genes_exons.table",sep="\t",header=T)
Klet_coverage<-read.table("Klet_arenosa_genes_exons.table",sep="\t",header=T)
Kowa_coverage<-read.table("Kowa_arenosa_genes_exons.table",sep="\t",header=T)
Wulm_coverage<-read.table("Wulm_arenosa_genes_exons.table",sep="\t",header=T)
Chok_coverage<-read.table("Chok_arenosa_genes_exons.table",sep="\t",header=T)

############################################################
#merge all populations and show all genes with coverage < 2#
############################################################
Coverage_all<-cbind(Mias_coverage,Zapa_coverage[,10:13],Klet_coverage[,10:13],Kowa_coverage[,10:13],Wulm_coverage[,10:13],Chok_coverage[,10:13])
names(Coverage_all)[10:33]<-c("Mias_mean_gene","Mias_median_gene","Mias_mean_exons","Mias_median_exons","Zapa_mean_gene","Zapa_median_gene","Zapa_mean_exons","Zapa_median_exons","Klet_mean_gene","Klet_median_gene","Klet_mean_exons","Klet_median_exons","Kowa_mean_gene","Kowa_median_gene","Kowa_mean_exons","Kowa_median_exons","Wulm_mean_gene","Wulm_median_gene","Wulm_mean_exons","Wulm_median_exons","Chok_mean_gene","Chok_median_gene","Chok_mean_exons","Chok_median_exons")
Numberofindividuals<-c(7,7,7,7,8,8,8,8,9,9,9,9,8,8,8,8,8,8,8,8,7,7,7,7)
Coverage_all[,10:33]<-Coverage_all[,10:33]/Numberofindividuals
Low_cov<-Coverage_all[apply(Coverage_all[,c(13,17,21,25,29,33)],MARGIN=1,function(x) all(x<4.7)),]
Low_cov_w_orthogroup= merge(Low_cov,Lyr_ortholist,by.x="ID",by.y="Gene",all.x=TRUE)
Low_cov_w_orthogroup_Thalgenes= merge(Low_cov_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Low_cov_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Low_cov_w_orthogroup_Thalgenes)[35]="Gene"
Thalgenes_described_Low_cov_w_orthogroup= merge(Thal_description, Low_cov_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Low_cov_all_genes<-merge(Thalgenes_described_Low_cov_w_orthogroup,Methomlist,by.x="Gene",by.y="AGI_code",all.x=TRUE)
write.table(Low_cov_all_genes,"Low_cov_metallist_arenosa_smaller4_7.table",sep="\t",row.names=F)
High_cov<-Coverage_all[apply(Coverage_all[,c(13,17,21,25,29,33)],MARGIN=1,function(x) all(x>44.8)),]
High_cov_w_orthogroup= merge(High_cov,Lyr_ortholist,by.x="ID",by.y="Gene",all.x=TRUE)
High_cov_w_orthogroup_Thalgenes= merge(High_cov_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(High_cov_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(High_cov_w_orthogroup_Thalgenes)[35]="Gene"
Thalgenes_described_High_cov_w_orthogroup= merge(Thal_description, High_cov_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
High_cov_all_genes<-merge(Thalgenes_described_High_cov_w_orthogroup,Methomlist,by.x="Gene",by.y="AGI_code",all.x=TRUE)
write.table(High_cov_all_genes,"High_cov_metallist_arenosa_bigger44_8.table",sep="\t",row.names=F)


#DP_MZ<-read.table("MiasZapaDP_woheader_arenosa.table")
#SNP DP: 9.4 - 22.4
############################################################
#GROM data#
############################################################
#lyrata

Mias_Indels<-read.table("Mias_Indels_arenosa.vcf",sep="\t",header=F)
Mias_CNVs1<-read.table("Mias_cnvs_arenosa.vcf",sep="\t",header=F)
Zapa_Indels<-read.table("Zapa_Indels_arenosa.vcf",sep="\t",header=F)
Zapa_CNVs1<-read.table("Zapa_cnvs_arenosa.vcf",sep="\t",header=F)
Klet_Indels<-read.table("Klet_Indels_arenosa.vcf",sep="\t",header=F)
Klet_CNVs1<-read.table("Klet_cnvs_arenosa.vcf",sep="\t",header=F)
Kowa_Indels<-read.table("Kowa_Indels_arenosa.vcf",sep="\t",header=F)
Kowa_CNVs1<-read.table("Kowa_cnvs_arenosa.vcf",sep="\t",header=F)
Wulm_Indels<-read.table("Wulm_Indels_arenosa.vcf",sep="\t",header=F)
Wulm_CNVs1<-read.table("Wulm_cnvs_arenosa.vcf",sep="\t",header=F)
Chok_Indels<-read.table("Chok_Indels_arenosa.vcf",sep="\t",header=F)
Chok_CNVs1<-read.table("Chok_cnvs_arenosa.vcf",sep="\t",header=F)

names(Mias_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Mias_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Zapa_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Zapa_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Klet_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Klet_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Kowa_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Kowa_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Wulm_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Wulm_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Chok_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Chok_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Scaffolds<-c("scaffold_1","scaffold_2","scaffold_3","scaffold_4","scaffold_5","scaffold_6","scaffold_7","scaffold_8")
Mias_CNVs1<-droplevels(Mias_CNVs1[Mias_CNVs1$Scaffold%in%Scaffolds,])
Zapa_CNVs1<-droplevels(Zapa_CNVs1[Zapa_CNVs1$Scaffold%in%Scaffolds,])
Klet_CNVs1<-droplevels(Klet_CNVs1[Klet_CNVs1$Scaffold%in%Scaffolds,])
Kowa_CNVs1<-droplevels(Kowa_CNVs1[Kowa_CNVs1$Scaffold%in%Scaffolds,])
Wulm_CNVs1<-droplevels(Wulm_CNVs1[Wulm_CNVs1$Scaffold%in%Scaffolds,])
Chok_CNVs1<-droplevels(Chok_CNVs1[Chok_CNVs1$Scaffold%in%Scaffolds,])

Mias_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Mias_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Mias_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Mias_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Mias_CNVs<-data.frame(Mias_CNVs1[,1:2],as.numeric(as.character(Mias_CNVs_END[,2])),Mias_CNVs1[,3:7],Mias_CNVs_info)
names(Mias_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Zapa_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Zapa_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Zapa_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Zapa_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Zapa_CNVs<-data.frame(Zapa_CNVs1[,1:2],as.numeric(as.character(Zapa_CNVs_END[,2])),Zapa_CNVs1[,3:7],Zapa_CNVs_info)
names(Zapa_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Klet_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Klet_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Klet_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Klet_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Klet_CNVs<-data.frame(Klet_CNVs1[,1:2],as.numeric(as.character(Klet_CNVs_END[,2])),Klet_CNVs1[,3:7],Klet_CNVs_info)
names(Klet_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Kowa_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Kowa_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Kowa_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Kowa_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Kowa_CNVs<-data.frame(Kowa_CNVs1[,1:2],as.numeric(as.character(Kowa_CNVs_END[,2])),Kowa_CNVs1[,3:7],Kowa_CNVs_info)
names(Kowa_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Wulm_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Wulm_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Wulm_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Wulm_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Wulm_CNVs<-data.frame(Wulm_CNVs1[,1:2],as.numeric(as.character(Wulm_CNVs_END[,2])),Wulm_CNVs1[,3:7],Wulm_CNVs_info)
names(Wulm_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Chok_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Chok_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Chok_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Chok_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Chok_CNVs<-data.frame(Chok_CNVs1[,1:2],as.numeric(as.character(Chok_CNVs_END[,2])),Chok_CNVs1[,3:7],Chok_CNVs_info)
names(Chok_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

#######################################################
#CN.mops data#
#######################################################
#change directory
Mias_cnvs_HC_strict<-read.csv("Mias_cnvs_only_arenosa.table",header=T,sep="\t")
Klet_cnvs_HC_strict<-read.csv("Klet_cnvs_only_arenosa.table",header=T,sep="\t")
Wulm_cnvs_HC_strict<-read.csv("Wulm_cnvs_only_arenosa.table",header=T,sep="\t")
Zapa_cnvs_HC_strict<-read.csv("Zapa_cnvs_only_arenosa.table",header=T,sep="\t")
Kowa_cnvs_HC_strict<-read.csv("Kowa_cnvs_only_arenosa.table",header=T,sep="\t")
Chok_cnvs_HC_strict<-read.csv("Chok_cnvs_only_arenosa.table",header=T,sep="\t")
#change directory

require(GenomicRanges)
#find overlapping regions between Mias and Zako with similar CNs

Mias_CNVs_GRange<-GRanges(seqnames=tolower(Mias_CNVs$Scaffold),ranges=IRanges(start=Mias_CNVs$Start_pos,end=Mias_CNVs$End_pos))
values(Mias_CNVs_GRange)<-Mias_CNVs[,4:12]
Zapa_CNVs_GRange<-GRanges(seqnames=tolower(Zapa_CNVs$Scaffold),ranges=IRanges(start=Zapa_CNVs$Start_pos,end=Zapa_CNVs$End_pos))
values(Zapa_CNVs_GRange)<-Zapa_CNVs[,4:12]
MZ_merge<-mergeByOverlaps(Mias_CNVs_GRange,Zapa_CNVs_GRange)
MZ_merge_same<-MZ_merge[abs(as.numeric(as.character(MZ_merge@ listData$ Mias_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(MZ_merge@ listData$ Zapa_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
MZ_merge_same_df<-as.data.frame(MZ_merge_same)[,1:3]
names(MZ_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
MZ_merge_same_df_GRange<-GRanges(seqnames=tolower(MZ_merge_same_df$Scaffold),ranges=IRanges(start=MZ_merge_same_df$Start_pos,end=MZ_merge_same_df$End_pos))

#subset Mias CNVs to regions without similar CN in Zapa

MZ_diff<-setdiff(Mias_CNVs_GRange,MZ_merge_same_df_GRange)
MZ<-subsetByOverlaps(Mias_CNVs_GRange,MZ_diff)
ZM_diff<-setdiff(Zapa_CNVs_GRange,MZ_merge_same_df_GRange)
ZM<-subsetByOverlaps(Zapa_CNVs_GRange,ZM_diff)

#overlap cn.mops and GROM regions

Mias_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Mias_cnvs_HC_strict$Chr),ranges=IRanges(start=Mias_cnvs_HC_strict$gene_start,end=Mias_cnvs_HC_strict$gene_end))
values(Mias_cnvs_HC_strict_GRange)<-cbind(Mias_cnvs_HC_strict[,1:7],Mias_cnvs_HC_strict[,11:37])
Mias_overlap_strict<-mergeByOverlaps(Mias_cnvs_HC_strict_GRange,MZ)
Mias_overlap_strict_df<-as.data.frame(Mias_overlap_strict)
Zapa_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Zapa_cnvs_HC_strict$Chr),ranges=IRanges(start=Zapa_cnvs_HC_strict$gene_start,end=Zapa_cnvs_HC_strict$gene_end))
values(Zapa_cnvs_HC_strict_GRange)<-cbind(Zapa_cnvs_HC_strict[,1:7],Zapa_cnvs_HC_strict[,11:37])
Zapa_overlap_strict<-mergeByOverlaps(Zapa_cnvs_HC_strict_GRange,ZM)
Zapa_overlap_strict_df<-as.data.frame(Zapa_overlap_strict)

Mias_CNVs_overlapping_strict<-data.frame(Mias_overlap_strict_df[,1:39],Mias_overlap_strict_df[,74:87])
names(Mias_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Gene","Type","Short_description","Curator_summary","Computational_description","OrthoGroup","Lyr_Gene","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene","gene_start","gene_end","gene_size","gene_strand","Widthofoverlap","ID","Name","Function","Localisation","EC.TC","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")
Zapa_CNVs_overlapping_strict<-data.frame(Zapa_overlap_strict_df[,1:39],Zapa_overlap_strict_df[,74:87])
names(Zapa_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Gene","Type","Short_description","Curator_summary","Computational_description","OrthoGroup","Lyr_Gene","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene","gene_start","gene_end","gene_size","gene_strand","Widthofoverlap","ID","Name","Function","Localisation","EC.TC","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Mias_CNVs_HC_GROM_strict<-rbind(Mias_CNVs_overlapping_strict,Zapa_CNVs_overlapping_strict)

CNVs_Mias_HC_GROM_cov<-merge(Mias_CNVs_HC_GROM_strict,Coverage_all,by.x="Lyr_Gene",by.y="ID",all.x=T)
CNVs_Mias_HC_GROM_cov<-CNVs_Mias_HC_GROM_cov[!((CNVs_Mias_HC_GROM_cov$CN_class=="CN0"|CNVs_Mias_HC_GROM_cov$CN_class=="CN1")&CNVs_Mias_HC_GROM_cov$GROM_ALT=="<DUP>"),]
CNVs_Mias_HC_GROM_cov<-CNVs_Mias_HC_GROM_cov[!(CNVs_Mias_HC_GROM_cov$GROM_ALT=="<DEL>"&!(CNVs_Mias_HC_GROM_cov$CN_class=="CN0"|CNVs_Mias_HC_GROM_cov$CN_class=="CN1")),]
write.table(CNVs_Mias_HC_GROM_cov,"Mias_CNVs_overlap_strict.table",row.names=F,sep="\t")


#find overlapping regions between Klet and Kowa with similar CNs

Klet_CNVs_GRange<-GRanges(seqnames=tolower(Klet_CNVs$Scaffold),ranges=IRanges(start=Klet_CNVs$Start_pos,end=Klet_CNVs$End_pos))
values(Klet_CNVs_GRange)<-Klet_CNVs[,4:12]
Kowa_CNVs_GRange<-GRanges(seqnames=tolower(Kowa_CNVs$Scaffold),ranges=IRanges(start=Kowa_CNVs$Start_pos,end=Kowa_CNVs$End_pos))
values(Kowa_CNVs_GRange)<-Kowa_CNVs[,4:12]
KK_merge<-mergeByOverlaps(Klet_CNVs_GRange,Kowa_CNVs_GRange)
KK_merge_same<-KK_merge[abs(as.numeric(as.character(KK_merge@ listData$ Klet_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(KK_merge@ listData$ Kowa_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
KK_merge_same_df<-as.data.frame(KK_merge_same)[,1:3]
names(KK_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
KK_merge_same_df_GRange<-GRanges(seqnames=tolower(KK_merge_same_df$Scaffold),ranges=IRanges(start=KK_merge_same_df$Start_pos,end=KK_merge_same_df$End_pos))

#subset Klet CNVs to regions without similar CN in Kowa

KK_diff<-setdiff(Klet_CNVs_GRange,KK_merge_same_df_GRange)
KK<-subsetByOverlaps(Klet_CNVs_GRange,KK_diff)
KoKl_diff<-setdiff(Kowa_CNVs_GRange,KK_merge_same_df_GRange)
KoKl<-subsetByOverlaps(Kowa_CNVs_GRange,KoKl_diff)

#overlap cn.mops and GROM regions

Klet_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Klet_cnvs_HC_strict$Chr),ranges=IRanges(start=Klet_cnvs_HC_strict$gene_start,end=Klet_cnvs_HC_strict$gene_end))
values(Klet_cnvs_HC_strict_GRange)<-cbind(Klet_cnvs_HC_strict[,1:7],Klet_cnvs_HC_strict[,11:37])
Klet_overlap_strict<-mergeByOverlaps(Klet_cnvs_HC_strict_GRange,KK)
Klet_overlap_strict_df<-as.data.frame(Klet_overlap_strict)
Kowa_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Kowa_cnvs_HC_strict$Chr),ranges=IRanges(start=Kowa_cnvs_HC_strict$gene_start,end=Kowa_cnvs_HC_strict$gene_end))
values(Kowa_cnvs_HC_strict_GRange)<-cbind(Kowa_cnvs_HC_strict[,1:7],Kowa_cnvs_HC_strict[,11:37])
Kowa_overlap_strict<-mergeByOverlaps(Kowa_cnvs_HC_strict_GRange,KoKl)
Kowa_overlap_strict_df<-as.data.frame(Kowa_overlap_strict)

Klet_CNVs_overlapping_strict<-data.frame(Klet_overlap_strict_df[,1:39],Klet_overlap_strict_df[,74:87])
names(Klet_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Gene","Type","Short_description","Curator_summary","Computational_description","OrthoGroup","Lyr_Gene","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene","gene_start","gene_end","gene_size","gene_strand","Widthofoverlap","ID","Name","Function","Localisation","EC.TC","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")
Kowa_CNVs_overlapping_strict<-data.frame(Kowa_overlap_strict_df[,1:39],Kowa_overlap_strict_df[,74:87])
names(Kowa_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Gene","Type","Short_description","Curator_summary","Computational_description","OrthoGroup","Lyr_Gene","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene","gene_start","gene_end","gene_size","gene_strand","Widthofoverlap","ID","Name","Function","Localisation","EC.TC","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Klet_CNVs_HC_GROM_strict<-rbind(Klet_CNVs_overlapping_strict,Kowa_CNVs_overlapping_strict)

CNVs_Klet_HC_GROM_cov<-merge(Klet_CNVs_HC_GROM_strict,Coverage_all,by.x="Lyr_Gene",by.y="ID",all.x=T)
CNVs_Klet_HC_GROM_cov<-CNVs_Klet_HC_GROM_cov[!((CNVs_Klet_HC_GROM_cov$CN_class=="CN0"|CNVs_Klet_HC_GROM_cov$CN_class=="CN1")&CNVs_Klet_HC_GROM_cov$GROM_ALT=="<DUP>"),]
CNVs_Klet_HC_GROM_cov<-CNVs_Klet_HC_GROM_cov[!(CNVs_Klet_HC_GROM_cov$GROM_ALT=="<DEL>"&!(CNVs_Klet_HC_GROM_cov$CN_class=="CN0"|CNVs_Klet_HC_GROM_cov$CN_class=="CN1")),]
write.table(CNVs_Klet_HC_GROM_cov,"Klet_CNVs_overlap_strict.table",row.names=F,sep="\t")

#find overlapping regions between Wulm and HD with similar CNs

Wulm_CNVs_GRange<-GRanges(seqnames=tolower(Wulm_CNVs$Scaffold),ranges=IRanges(start=Wulm_CNVs$Start_pos,end=Wulm_CNVs$End_pos))
values(Wulm_CNVs_GRange)<-Wulm_CNVs[,4:12]
Chok_CNVs_GRange<-GRanges(seqnames=tolower(Chok_CNVs$Scaffold),ranges=IRanges(start=Chok_CNVs$Start_pos,end=Chok_CNVs$End_pos))
values(Chok_CNVs_GRange)<-Chok_CNVs[,4:12]
WC_merge<-mergeByOverlaps(Wulm_CNVs_GRange,Chok_CNVs_GRange)
WC_merge_same<-WC_merge[abs(as.numeric(as.character(WC_merge@ listData$ Wulm_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(WC_merge@ listData$ Chok_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
WC_merge_same_df<-as.data.frame(WC_merge_same)[,1:3]
names(WC_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
WC_merge_same_df_GRange<-GRanges(seqnames=tolower(WC_merge_same_df$Scaffold),ranges=IRanges(start=WC_merge_same_df$Start_pos,end=WC_merge_same_df$End_pos))

#subset Wulm CNVs to regions without similar CN in Chok

WC_diff<-setdiff(Wulm_CNVs_GRange,WC_merge_same_df_GRange)
WC<-subsetByOverlaps(Wulm_CNVs_GRange,WC_diff)
CW_diff<-setdiff(Chok_CNVs_GRange,WC_merge_same_df_GRange)
CW<-subsetByOverlaps(Chok_CNVs_GRange,CW_diff)

#overlap cn.mops and GROM regions

Wulm_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Wulm_cnvs_HC_strict$Chr),ranges=IRanges(start=Wulm_cnvs_HC_strict$gene_start,end=Wulm_cnvs_HC_strict$gene_end))
values(Wulm_cnvs_HC_strict_GRange)<-cbind(Wulm_cnvs_HC_strict[,1:7],Wulm_cnvs_HC_strict[,11:37])
Wulm_overlap_strict<-mergeByOverlaps(Wulm_cnvs_HC_strict_GRange,WC)
Wulm_overlap_strict_df<-as.data.frame(Wulm_overlap_strict)
Chok_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Chok_cnvs_HC_strict$Chr),ranges=IRanges(start=Chok_cnvs_HC_strict$gene_start,end=Chok_cnvs_HC_strict$gene_end))
values(Chok_cnvs_HC_strict_GRange)<-cbind(Chok_cnvs_HC_strict[,1:7],Chok_cnvs_HC_strict[,11:37])
Chok_overlap_strict<-mergeByOverlaps(Chok_cnvs_HC_strict_GRange,CW)
Chok_overlap_strict_df<-as.data.frame(Chok_overlap_strict)

Wulm_CNVs_overlapping_strict<-data.frame(Wulm_overlap_strict_df[,1:39],Wulm_overlap_strict_df[,74:87])
names(Wulm_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Gene","Type","Short_description","Curator_summary","Computational_description","OrthoGroup","Lyr_Gene","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene","gene_start","gene_end","gene_size","gene_strand","Widthofoverlap","ID","Name","Function","Localisation","EC.TC","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")
Chok_CNVs_overlapping_strict<-data.frame(Chok_overlap_strict_df[,1:39],Chok_overlap_strict_df[,74:87])
names(Chok_CNVs_overlapping_strict)<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand","Gene","Type","Short_description","Curator_summary","Computational_description","OrthoGroup","Lyr_Gene","Copy_width","Copy_strand","sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene","gene_start","gene_end","gene_size","gene_strand","Widthofoverlap","ID","Name","Function","Localisation","EC.TC","Comment","BINCODE","BIN","EvidenceCode","Annotator...Curator","metal.1","metal.2","metal.3",
"GROM_scaffold","GROM_Start_pos","GROM_End_pos","GROM_width","GROM_strand","GROM_ID","GROM_REF","GROM_ALT","GROM_QUAL","GROM_FILTER","CNV_standard_deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Wulm_CNVs_HC_GROM_strict<-rbind(Wulm_CNVs_overlapping_strict,Chok_CNVs_overlapping_strict)

CNVs_Wulm_HC_GROM_cov<-merge(Wulm_CNVs_HC_GROM_strict,Coverage_all,by.x="Lyr_Gene",by.y="ID",all.x=T)
CNVs_Wulm_HC_GROM_cov<-CNVs_Wulm_HC_GROM_cov[!((CNVs_Wulm_HC_GROM_cov$CN_class=="CN0"|CNVs_Wulm_HC_GROM_cov$CN_class=="CN1")&CNVs_Wulm_HC_GROM_cov$GROM_ALT=="<DUP>"),]
CNVs_Wulm_HC_GROM_cov<-CNVs_Wulm_HC_GROM_cov[!(CNVs_Wulm_HC_GROM_cov$GROM_ALT=="<DEL>"&!(CNVs_Wulm_HC_GROM_cov$CN_class=="CN0"|CNVs_Wulm_HC_GROM_cov$CN_class=="CN1")),]
write.table(CNVs_Wulm_HC_GROM_cov,"Wulm_CNVs_overlap_strict.table",row.names=F,sep="\t")

###################################################
#INDEL analysis#
###################################################
require(xlsx)

############################################################
#GROM data#
############################################################
#lyrata

Klet_indel_grom<-read.table("KletKowa_INDEL_GROM.arenosa.ann.vcf",sep="\t",header=F)
Kowa_indel_grom<-read.table("KletKowa_INDEL_GROM.Kowa.arenosa.ann.vcf",sep="\t",header=F)
names(Klet_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Kowa_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Mias_indel_grom<-read.table("MiasZapa_INDEL_GROM.arenosa.ann.vcf",sep="\t",header=F)
Zapa_indel_grom<-read.table("MiasZapa_INDEL_GROM.Zapa.arenosa.ann.vcf",sep="\t",header=F)
names(Mias_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Zapa_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

Wulm_indel_grom<-read.table("WulmHD_INDEL_GROM.ann.vcf",sep="\t",header=F)
HD_indel_grom<-read.table("WulmHD_INDEL_GROM.HD.ann.vcf",sep="\t",header=F)
names(Wulm_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(HD_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")


#change workingdirectory
MZ_INDEL<-read.table("SNP_MZ_Indels_diffabove04_hQD_arenosa.txt",sep="\t",header=T)
MZ_INDEL[,8]<-as.character(MZ_INDEL[,8])
MZ_INDEL[,9]<-as.numeric(as.character(MZ_INDEL[,9]))

testM_Indels<-merge(MZ_INDEL,Mias_Indels,by.x=c("Chrom","Pos"),by.y=c("Scaffold","Pos"))
write.table(testM_Indels,"AllMias_Indels_overlap.table",row.names=F,sep="\t")

KK_INDEL<-read.table("SNP_KK_Indels_diffabove04_hQD_arenosa.txt",sep="\t",header=T)
KK_INDEL[,8]<-as.character(KK_INDEL[,8])
KK_INDEL[,9]<-as.numeric(as.character(KK_INDEL[,9]))

KletKowa_Indels<-rbind(Klet_indel_grom,Kowa_indel_grom)
KK_INDEL_GR<-GRanges(seqnames=tolower(KK_INDEL$Chrom),ranges=IRanges(start=KK_INDEL$Pos,end=KK_INDEL$Pos))
KletKowa_Indels_GR<-GRanges(seqnames=tolower(KletKowa_Indels$Scaffold),ranges=IRanges(start=KletKowa_Indels$Pos,end=KletKowa_Indels$Pos))
values(KK_INDEL_GR)<-KK_INDEL[,c(1:7,10:29)]
values(KletKowa_Indels_GR)<-KletKowa_Indels[,3:10]
KK_Indels_nearest<-nearest(KK_INDEL_GR,KletKowa_Indels_GR,ignore.strand=T)
KK_nearest<-KletKowa_Indels[KK_Indels_nearest,1:10]
KK_Indel_table<-cbind(KK_INDEL,KK_nearest)
KK_Indel_table2<-KK_Indel_table[abs(KK_Indel_table[,9]-KK_Indel_table[,31])<=max(nchar(as.character(KK_Indel_table[,12])),nchar(as.character(KK_Indel_table[,13])),nchar(as.character(KK_Indel_table[,33])),nchar(as.character(KK_Indel_table[,34]))),]

SNPeffAnnsplita<-strsplit(sub(";","bla",as.character(KK_Indel_table2$INFO)),"bla",fixed=T)
SNPeffAnnsplitb<-lapply(SNPeffAnnsplita,"[",2:max(unlist(lapply(SNPeffAnnsplita,length))))
SNPeffAnnsplit<-strsplit(as.character(SNPeffAnnsplitb),",",fixed=T)
SNPeffAnnsplit1<-lapply(SNPeffAnnsplit,"[",1:max(unlist(lapply(SNPeffAnnsplit,length))))
SNPeffAnnsplit2<-data.frame(matrix(unlist(SNPeffAnnsplit1),nrow=length(SNPeffAnnsplit1),byrow=T))
SNPeffdf<-data.frame()
for (i in 1:max(unlist(lapply(SNPeffAnnsplit,length))))
	{SNPeff1<-ifelse(!is.na(SNPeffAnnsplit2[,i]),strsplit(as.character(SNPeffAnnsplit2[,i]),"|",fixed=T),NA)
	SNPeffAnn1<-lapply(SNPeff1,"[",1:4)
	SNPeffAnndf1<-data.frame(matrix(unlist(SNPeffAnn1),nrow=length(SNPeffAnnsplit1),byrow=T))
	names(SNPeffAnndf1)<-c("Base","Characterization","Effect","Gene")
	SNPeffAnndf2<-ifelse(!is.na(SNPeffAnndf1[,1]),strsplit(as.character(SNPeffAnndf1[,1]),"=",fixed=T),NA)
	SNPeffAnndf3<-lapply(SNPeffAnndf2,"[",2:ifelse(max(unlist(lapply(SNPeffAnndf2,length)))==1,2,max(unlist(lapply(SNPeffAnndf2,length)))))
	SNPeffAnndf1$Base<-c(unlist(SNPeffAnndf3))
	SNPeffcomp1<-cbind(KK_Indel_table2[,c(1:36,38:39)],SNPeffAnndf1)
	SNPeffdf<-rbind(SNPeffdf,SNPeffcomp1[!is.na(SNPeffcomp1[,39]),])
	}
SNPeffdf<-SNPeffdf[,1:41]

SNPeffdf2<-SNPeffdf[SNPeffdf[,41]=="HIGH"&SNPeffdf[,16]=="HIGH",]
SNPeffdf3<-SNPeffdf[SNPeffdf[,41]=="HIGH"|SNPeffdf[,16]=="HIGH",]
nrow(SNPeffdf2[!duplicated(SNPeffdf2),])
nrow(SNPeffdf3[!duplicated(SNPeffdf3),])
SNPeffdf2<-SNPeffdf2[!duplicated(SNPeffdf2),]
SNPeffdf3<-SNPeffdf3[!duplicated(SNPeffdf3),]

#remove duplicates where different high effect mutations occurred in both populations
dup_genes<-SNPeffdf2$Lyr_Gene[duplicated(SNPeffdf2$Lyr_Gene)]
unique_pos<-!(duplicated(SNPeffdf2$Pos))
Gene_tab<-tabulate(SNPeffdf2$Lyr_Gene)
sum_AF<-sum(SNPeffdf2$AF_Klet>SNPeffdf2$AF_Kowa)

removal_df<-SNPeffdf2[SNPeffdf2$Lyr_Gene%in%dup_genes&unique_pos=="TRUE",]
removal<-c()
for(i in removal_df$Lyr_Gene)
	{for (j in unique(removal_df$Gene[removal_df$Lyr_Gene==i]))
		{if (is.na(j)==TRUE)
			{removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i,])!=sum(SNPeffdf2$AF_Klet[SNPeffdf2$Lyr_Gene==i]>SNPeffdf2$AF_Kowa[SNPeffdf2$Lyr_Gene==i])&nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i,])!=sum(SNPeffdf2$AF_Klet[SNPeffdf2$Lyr_Gene==i]<SNPeffdf2$AF_Kowa[SNPeffdf2$Lyr_Gene==i]),paste0(i),"keep"))
			}
		else
			{
			removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j,])!=sum(SNPeffdf2$AF_Klet[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j]>SNPeffdf2$AF_Kowa[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j])&nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j,])!=sum(SNPeffdf2$AF_Klet[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j]<SNPeffdf2$AF_Kowa[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j]),paste0(i,"_",j),"keep"))
			}
		}
	}

SNPeffdf2_clean<-SNPeffdf2[!(paste0(SNPeffdf2$Lyr_Gene,"_",SNPeffdf2$Gene)%in%removal|(is.na(SNPeffdf2$Gene)&SNPeffdf2$Lyr_Gene%in%removal)),]

SNPeffdf2_clean<-SNPeffdf2_clean[!duplicated(SNPeffdf2_clean[,c(1,6,7,8,9,12,13,30,31,33,34)]),]
write.table(SNPeffdf2_clean,"Indel_KK_HC_GROM_arenosa.table",row.names=F,sep="\t")

MZ_INDEL<-read.table("SNP_MZ_Indels_diffabove04_hQD_arenosa.txt",sep="\t",header=T)
MZ_INDEL[,8]<-as.character(MZ_INDEL[,8])
MZ_INDEL[,9]<-as.numeric(as.character(MZ_INDEL[,9]))

MiasZapa_Indels<-rbind(Mias_indel_grom,Zapa_indel_grom)
MZ_INDEL_GR<-GRanges(seqnames=tolower(MZ_INDEL$Chrom),ranges=IRanges(start=MZ_INDEL$Pos,end=MZ_INDEL$Pos))
MiasZapa_Indels_GR<-GRanges(seqnames=tolower(MiasZapa_Indels$Scaffold),ranges=IRanges(start=MiasZapa_Indels$Pos,end=MiasZapa_Indels$Pos))
values(MZ_INDEL_GR)<-MZ_INDEL[,c(1:7,10:29)]
values(MiasZapa_Indels_GR)<-MiasZapa_Indels[,3:10]
MZ_Indels_nearest<-nearest(MZ_INDEL_GR,MiasZapa_Indels_GR,ignore.strand=T)
MZ_nearest<-MiasZapa_Indels[MZ_Indels_nearest,1:10]
MZ_Indel_table<-cbind(MZ_INDEL,MZ_nearest)
MZ_Indel_table2<-MZ_Indel_table[abs(MZ_Indel_table[,9]-MZ_Indel_table[,31])<=max(nchar(as.character(MZ_Indel_table[,12])),nchar(as.character(MZ_Indel_table[,13])),nchar(as.character(MZ_Indel_table[,33])),nchar(as.character(MZ_Indel_table[,34]))),]

SNPeffAnnsplita<-strsplit(sub(";","bla",as.character(MZ_Indel_table2$INFO)),"bla",fixed=T)
SNPeffAnnsplitb<-lapply(SNPeffAnnsplita,"[",2:max(unlist(lapply(SNPeffAnnsplita,length))))
SNPeffAnnsplit<-strsplit(as.character(SNPeffAnnsplitb),",",fixed=T)
SNPeffAnnsplit1<-lapply(SNPeffAnnsplit,"[",1:max(unlist(lapply(SNPeffAnnsplit,length))))
SNPeffAnnsplit2<-data.frame(matrix(unlist(SNPeffAnnsplit1),nrow=length(SNPeffAnnsplit1),byrow=T))
SNPeffdf<-data.frame()
for (i in 1:max(unlist(lapply(SNPeffAnnsplit,length))))
	{SNPeff1<-ifelse(!is.na(SNPeffAnnsplit2[,i]),strsplit(as.character(SNPeffAnnsplit2[,i]),"|",fixed=T),NA)
	SNPeffAnn1<-lapply(SNPeff1,"[",1:4)
	SNPeffAnndf1<-data.frame(matrix(unlist(SNPeffAnn1),nrow=length(SNPeffAnnsplit1),byrow=T))
	names(SNPeffAnndf1)<-c("Base","Characterization","Effect","Gene")
	SNPeffAnndf2<-ifelse(!is.na(SNPeffAnndf1[,1]),strsplit(as.character(SNPeffAnndf1[,1]),"=",fixed=T),NA)
	SNPeffAnndf3<-lapply(SNPeffAnndf2,"[",2:ifelse(max(unlist(lapply(SNPeffAnndf2,length)))==1,2,max(unlist(lapply(SNPeffAnndf2,length)))))
	SNPeffAnndf1$Base<-c(unlist(SNPeffAnndf3))
	SNPeffcomp1<-cbind(MZ_Indel_table2[,c(1:36,38:39)],SNPeffAnndf1)
	SNPeffdf<-rbind(SNPeffdf,SNPeffcomp1[!is.na(SNPeffcomp1[,39]),])
	}
SNPeffdf<-SNPeffdf[,1:41]

SNPeffdf2<-SNPeffdf[SNPeffdf[,41]=="HIGH"&SNPeffdf[,16]=="HIGH",]
SNPeffdf3<-SNPeffdf[SNPeffdf[,41]=="HIGH"|SNPeffdf[,16]=="HIGH",]
nrow(SNPeffdf2[!duplicated(SNPeffdf2),])
nrow(SNPeffdf3[!duplicated(SNPeffdf3),])
SNPeffdf2<-SNPeffdf2[!duplicated(SNPeffdf2),]
SNPeffdf3<-SNPeffdf3[!duplicated(SNPeffdf3),]

#remove duplicates where different high effect mutations occurred in both populations
dup_genes<-SNPeffdf2$Lyr_Gene[duplicated(SNPeffdf2$Lyr_Gene)]
unique_pos<-!(duplicated(SNPeffdf2$Pos))
Gene_tab<-tabulate(SNPeffdf2$Lyr_Gene)
sum_AF<-sum(SNPeffdf2$AF_Mias>SNPeffdf2$AF_Zako)

removal_df<-SNPeffdf2[SNPeffdf2$Lyr_Gene%in%dup_genes&unique_pos=="TRUE",]
removal<-c()
for(i in removal_df$Lyr_Gene)
	{for (j in unique(removal_df$Gene[removal_df$Lyr_Gene==i]))
		{if (is.na(j)==TRUE)
			{removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i,])!=sum(SNPeffdf2$AF_Mias[SNPeffdf2$Lyr_Gene==i]>SNPeffdf2$AF_Zako[SNPeffdf2$Lyr_Gene==i])&nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i,])!=sum(SNPeffdf2$AF_Mias[SNPeffdf2$Lyr_Gene==i]<SNPeffdf2$AF_Zako[SNPeffdf2$Lyr_Gene==i]),paste0(i),"keep"))
			}
		else
			{
			removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j,])!=sum(SNPeffdf2$AF_Mias[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j]>SNPeffdf2$AF_Zako[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j])&nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j,])!=sum(SNPeffdf2$AF_Mias[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j]<SNPeffdf2$AF_Zako[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j]),paste0(i,"_",j),"keep"))
			}
		}
	}

SNPeffdf2_clean<-SNPeffdf2[!(paste0(SNPeffdf2$Lyr_Gene,"_",SNPeffdf2$Gene)%in%removal|(is.na(SNPeffdf2$Gene)&SNPeffdf2$Lyr_Gene%in%removal)),]

SNPeffdf2_clean<-SNPeffdf2_clean[!duplicated(SNPeffdf2_clean[,c(1,6,7,8,9,12,13,30,31,33,34)]),]
write.table(SNPeffdf2_clean,"Indel_MZ_HC_GROM_arenosa.table",row.names=F,sep="\t")



WH_INDEL<-read.table("SNP_WH_Indels_diffabove04_hQD.txt",sep="\t",header=T)
WH_INDEL[,8]<-as.character(WH_INDEL[,8])
WH_INDEL[,9]<-as.numeric(as.character(WH_INDEL[,9]))

WulmHD_Indels<-rbind(Wulm_indel_grom,HD_indel_grom)
WH_INDEL_GR<-GRanges(seqnames=tolower(WH_INDEL$Chrom),ranges=IRanges(start=WH_INDEL$Pos,end=WH_INDEL$Pos))
WulmHD_Indels_GR<-GRanges(seqnames=tolower(WulmHD_Indels$Scaffold),ranges=IRanges(start=WulmHD_Indels$Pos,end=WulmHD_Indels$Pos))
values(WH_INDEL_GR)<-WH_INDEL[,c(1:7,10:29)]
values(WulmHD_Indels_GR)<-WulmHD_Indels[,3:10]
WH_Indels_nearest<-nearest(WH_INDEL_GR,WulmHD_Indels_GR,ignore.strand=T)
WH_nearest<-WulmHD_Indels[WH_Indels_nearest,1:10]
WH_Indel_table<-cbind(WH_INDEL,WH_nearest)
WH_Indel_table2<-WH_Indel_table[abs(WH_Indel_table[,9]-WH_Indel_table[,31])<=max(nchar(as.character(WH_Indel_table[,12])),nchar(as.character(WH_Indel_table[,13])),nchar(as.character(WH_Indel_table[,33])),nchar(as.character(WH_Indel_table[,34]))),]

SNPeffAnnsplita<-strsplit(sub(";","bla",as.character(WH_Indel_table2$INFO)),"bla",fixed=T)
SNPeffAnnsplitb<-lapply(SNPeffAnnsplita,"[",2:max(unlist(lapply(SNPeffAnnsplita,length))))
SNPeffAnnsplit<-strsplit(as.character(SNPeffAnnsplitb),",",fixed=T)
SNPeffAnnsplit1<-lapply(SNPeffAnnsplit,"[",1:max(unlist(lapply(SNPeffAnnsplit,length))))
SNPeffAnnsplit2<-data.frame(matrix(unlist(SNPeffAnnsplit1),nrow=length(SNPeffAnnsplit1),byrow=T))
SNPeffdf<-data.frame()
for (i in 1:max(unlist(lapply(SNPeffAnnsplit,length))))
	{SNPeff1<-ifelse(!is.na(SNPeffAnnsplit2[,i]),strsplit(as.character(SNPeffAnnsplit2[,i]),"|",fixed=T),NA)
	SNPeffAnn1<-lapply(SNPeff1,"[",1:4)
	SNPeffAnndf1<-data.frame(matrix(unlist(SNPeffAnn1),nrow=length(SNPeffAnnsplit1),byrow=T))
	names(SNPeffAnndf1)<-c("Base","Characterization","Effect","Gene")
	SNPeffAnndf2<-ifelse(!is.na(SNPeffAnndf1[,1]),strsplit(as.character(SNPeffAnndf1[,1]),"=",fixed=T),NA)
	SNPeffAnndf3<-lapply(SNPeffAnndf2,"[",2:ifelse(max(unlist(lapply(SNPeffAnndf2,length)))==1,2,max(unlist(lapply(SNPeffAnndf2,length)))))
	SNPeffAnndf1$Base<-c(unlist(SNPeffAnndf3))
	SNPeffcomp1<-cbind(WH_Indel_table2[,c(1:36,38:39)],SNPeffAnndf1)
	SNPeffdf<-rbind(SNPeffdf,SNPeffcomp1[!is.na(SNPeffcomp1[,39]),])
	}
SNPeffdf<-SNPeffdf[,1:41]

SNPeffdf2<-SNPeffdf[SNPeffdf[,41]=="HIGH"&SNPeffdf[,16]=="HIGH",]
SNPeffdf3<-SNPeffdf[SNPeffdf[,41]=="HIGH"|SNPeffdf[,16]=="HIGH",]
nrow(SNPeffdf2[!duplicated(SNPeffdf2),])
nrow(SNPeffdf3[!duplicated(SNPeffdf3),])
SNPeffdf2<-SNPeffdf2[!duplicated(SNPeffdf2),]
SNPeffdf3<-SNPeffdf3[!duplicated(SNPeffdf3),]

#remove duplicates where different high effect mutations occurred in both populations
dup_genes<-SNPeffdf2$Lyr_Gene[duplicated(SNPeffdf2$Lyr_Gene)]
unique_pos<-!(duplicated(SNPeffdf2$Pos))
Gene_tab<-tabulate(SNPeffdf2$Lyr_Gene)
sum_AF<-sum(SNPeffdf2$AF_Wulm>SNPeffdf2$AF_HD)

removal_df<-SNPeffdf2[SNPeffdf2$Lyr_Gene%in%dup_genes&unique_pos=="TRUE",]
removal<-c()
for(i in removal_df$Lyr_Gene)
	{for (j in unique(removal_df$Gene[removal_df$Lyr_Gene==i]))
		{if (is.na(j)==TRUE)
			{removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i,])!=sum(SNPeffdf2$AF_Wulm[SNPeffdf2$Lyr_Gene==i]>SNPeffdf2$AF_HD[SNPeffdf2$Lyr_Gene==i])&nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i,])!=sum(SNPeffdf2$AF_Wulm[SNPeffdf2$Lyr_Gene==i]<SNPeffdf2$AF_HD[SNPeffdf2$Lyr_Gene==i]),paste0(i),"keep"))
			}
		else
			{
			removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j,])!=sum(SNPeffdf2$AF_Wulm[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j]>SNPeffdf2$AF_HD[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j])&nrow(SNPeffdf2[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j,])!=sum(SNPeffdf2$AF_Wulm[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j]<SNPeffdf2$AF_HD[SNPeffdf2$Lyr_Gene==i&SNPeffdf2$Gene==j]),paste0(i,"_",j),"keep"))
			}
		}
	}

SNPeffdf2_clean<-SNPeffdf2[!(paste0(SNPeffdf2$Lyr_Gene,"_",SNPeffdf2$Gene)%in%removal|(is.na(SNPeffdf2$Gene)&SNPeffdf2$Lyr_Gene%in%removal)),]

SNPeffdf2_clean<-SNPeffdf2_clean[!duplicated(SNPeffdf2_clean[,c(1,6,7,8,9,12,13,30,31,33,34)]),]
write.table(SNPeffdf2_clean,"Indel_WH_HC_GROM.table",row.names=F,sep="\t")







