require(xlsx)

FST<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",1)
DXY<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",2)
AFDabs<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",3)
NIELSEN<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",4)
VARLD<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",5)
FLK<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",6)
DD<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",7)

TajimasD_Klet<-read.table("TajimasDKletlyratamean_arenosa.csv",header=FALSE,sep="\t")
names(TajimasD_Klet)<-c("Scaffold","Start_pos","End_pos","TajimasDKlet")

TajimasD_Kowa<-read.table("TajimasDKowalyratamean_arenosa.csv",header=FALSE,sep="\t")
names(TajimasD_Kowa)<-c("Scaffold","Start_pos","End_pos","TajimasDKowa")

FayandWusH_Klet<-read.table("FayandWusHKletlyrata_mean_arenosa.csv",header=FALSE,sep="\t")
names(FayandWusH_Klet)<-c("Scaffold","Start_pos","End_pos","FayandWusHKlet")

FayandWusH_Kowa<-read.table("FayandWusHKowalyrata_mean_arenosa.csv",header=FALSE,sep="\t")
names(FayandWusH_Kowa)<-c("Scaffold","Start_pos","End_pos","FayandWusHKowa")

FST$Window_start<-as.numeric(as.character(FST$Window_start))
FST$Window_end<-as.numeric(as.character(FST$Window_end))
DXY$Window_start<-as.numeric(as.character(DXY$Window_start))
DXY$Window_end<-as.numeric(as.character(DXY$Window_end))
AFDabs$Window_start<-as.numeric(as.character(AFDabs$Window_start))
AFDabs$Window_end<-as.numeric(as.character(AFDabs$Window_end))
NIELSEN$Window_start<-as.numeric(as.character(NIELSEN$Window_start))
NIELSEN$Window_end<-as.numeric(as.character(NIELSEN$Window_end))
VARLD$Window_start<-as.numeric(as.character(VARLD$Window_start))
VARLD$Window_end<-as.numeric(as.character(VARLD$Window_end))
FLK$Window_start<-as.numeric(as.character(FLK$Window_start))
FLK$Window_end<-as.numeric(as.character(FLK$Window_end))
DD$Window_start<-as.numeric(as.character(DD$Window_start))
DD$Window_end<-as.numeric(as.character(DD$Window_end))

require(GenomicRanges)
FST_GRange<-GRanges(seqnames=tolower(FST$Chr),ranges=IRanges(start=FST$Window_start,end=FST$Window_end))
DXY_GRange<-GRanges(seqnames=tolower(DXY$Chr),ranges=IRanges(start=DXY$Window_start,end=DXY$Window_end))
AFDabs_GRange<-GRanges(seqnames=tolower(AFDabs$Chr),ranges=IRanges(start=AFDabs$Window_start,end=AFDabs$Window_end))
NIELSEN_GRange<-GRanges(seqnames=tolower(NIELSEN$Chr),ranges=IRanges(start=NIELSEN$Window_start,end=NIELSEN$Window_end))
VARLD_GRange<-GRanges(seqnames=tolower(VARLD$Chr),ranges=IRanges(start=VARLD$Window_start,end=VARLD$Window_end))
FLK_GRange<-GRanges(seqnames=tolower(FLK$Chr),ranges=IRanges(start=FLK$Window_start,end=FLK$Window_end))
DD_GRange<-GRanges(seqnames=tolower(DD$Chr),ranges=IRanges(start=DD$Window_start,end=DD$Window_end))

values(FST_GRange)<-data.frame(FST[,1:7],FST[,11:36])
values(DXY_GRange)<-data.frame(DXY[,1:7],DXY[,11:36])
values(AFDabs_GRange)<-data.frame(AFDabs[,1:7],AFDabs[,11:36])
values(NIELSEN_GRange)<-data.frame(NIELSEN[,1:7],NIELSEN[,11:36])
values(VARLD_GRange)<-data.frame(VARLD[,1:7],VARLD[,11:36])
values(FLK_GRange)<-data.frame(FLK[,1:7],FLK[,11:36])
values(DD_GRange)<-data.frame(DD[,1:7],DD[,11:36])

TajimasD_Klet_GRange<-GRanges(seqnames=tolower(TajimasD_Klet$Scaffold),ranges=IRanges(start=TajimasD_Klet$Start_pos,end=TajimasD_Klet$End_pos))
TajimasD_Kowa_GRange<-GRanges(seqnames=tolower(TajimasD_Kowa$Scaffold),ranges=IRanges(start=TajimasD_Kowa$Start_pos,end=TajimasD_Kowa$End_pos))
FayandWusH_Klet_GRange<-GRanges(seqnames=tolower(FayandWusH_Klet$Scaffold),ranges=IRanges(start=FayandWusH_Klet$Start_pos,end=FayandWusH_Klet$End_pos))
FayandWusH_Kowa_GRange<-GRanges(seqnames=tolower(FayandWusH_Kowa$Scaffold),ranges=IRanges(start=FayandWusH_Kowa$Start_pos,end=FayandWusH_Kowa$End_pos))

values(TajimasD_Klet_GRange)<-TajimasD_Klet[,4]
values(TajimasD_Kowa_GRange)<-TajimasD_Kowa[,4]
values(FayandWusH_Klet_GRange)<-FayandWusH_Klet[,4]
values(FayandWusH_Kowa_GRange)<-FayandWusH_Kowa[,4]

TD_Klet_nearestFST<-nearest(FST_GRange,TajimasD_Klet_GRange,ignore.strand=T)
FST$TD_Klet_nearest<-TajimasD_Klet[TD_Klet_nearestFST,4]
TD_Klet_nearestDXY<-nearest(DXY_GRange,TajimasD_Klet_GRange,ignore.strand=T)
DXY$TD_Klet_nearest<-TajimasD_Klet[TD_Klet_nearestDXY,4]
TD_Klet_nearestAFDabs<-nearest(AFDabs_GRange,TajimasD_Klet_GRange,ignore.strand=T)
AFDabs$TD_Klet_nearest<-TajimasD_Klet[TD_Klet_nearestAFDabs,4]
TD_Klet_nearestNIELSEN<-nearest(NIELSEN_GRange,TajimasD_Klet_GRange,ignore.strand=T)
NIELSEN$TD_Klet_nearest<-TajimasD_Klet[TD_Klet_nearestNIELSEN,4]
TD_Klet_nearestVARLD<-nearest(VARLD_GRange,TajimasD_Klet_GRange,ignore.strand=T)
VARLD$TD_Klet_nearest<-TajimasD_Klet[TD_Klet_nearestVARLD,4]
TD_Klet_nearestFLK<-nearest(FLK_GRange,TajimasD_Klet_GRange,ignore.strand=T)
FLK$TD_Klet_nearest<-TajimasD_Klet[TD_Klet_nearestFLK,4]
TD_Klet_nearestDD<-nearest(DD_GRange,TajimasD_Klet_GRange,ignore.strand=T)
DD$TD_Klet_nearest<-TajimasD_Klet[TD_Klet_nearestDD,4]


#test<-data.frame(TajimasD_Klet[TD_Klet_nearest,],FST)
#test2<-data.frame(test$Chr,test$Window_start,test$Window_end,test$Scaffold,test$Start_pos,test$End_pos)

TD_Kowa_nearestFST<-nearest(FST_GRange,TajimasD_Kowa_GRange,ignore.strand=T)
FST$TD_Kowa_nearest<-TajimasD_Kowa[TD_Kowa_nearestFST,4]
TD_Kowa_nearestDXY<-nearest(DXY_GRange,TajimasD_Kowa_GRange,ignore.strand=T)
DXY$TD_Kowa_nearest<-TajimasD_Kowa[TD_Kowa_nearestDXY,4]
TD_Kowa_nearestAFDabs<-nearest(AFDabs_GRange,TajimasD_Kowa_GRange,ignore.strand=T)
AFDabs$TD_Kowa_nearest<-TajimasD_Kowa[TD_Kowa_nearestAFDabs,4]
TD_Kowa_nearestNIELSEN<-nearest(NIELSEN_GRange,TajimasD_Kowa_GRange,ignore.strand=T)
NIELSEN$TD_Kowa_nearest<-TajimasD_Kowa[TD_Kowa_nearestNIELSEN,4]
TD_Kowa_nearestVARLD<-nearest(VARLD_GRange,TajimasD_Kowa_GRange,ignore.strand=T)
VARLD$TD_Kowa_nearest<-TajimasD_Kowa[TD_Kowa_nearestVARLD,4]
TD_Kowa_nearestFLK<-nearest(FLK_GRange,TajimasD_Kowa_GRange,ignore.strand=T)
FLK$TD_Kowa_nearest<-TajimasD_Kowa[TD_Kowa_nearestFLK,4]
TD_Kowa_nearestDD<-nearest(DD_GRange,TajimasD_Kowa_GRange,ignore.strand=T)
DD$TD_Kowa_nearest<-TajimasD_Kowa[TD_Kowa_nearestDD,4]

FWH_Klet_nearestFST<-nearest(FST_GRange,FayandWusH_Klet_GRange,ignore.strand=T)
FST$FWH_Klet_nearest<-FayandWusH_Klet[FWH_Klet_nearestFST,4]
FWH_Klet_nearestDXY<-nearest(DXY_GRange,FayandWusH_Klet_GRange,ignore.strand=T)
DXY$FWH_Klet_nearest<-FayandWusH_Klet[FWH_Klet_nearestDXY,4]
FWH_Klet_nearestAFDabs<-nearest(AFDabs_GRange,FayandWusH_Klet_GRange,ignore.strand=T)
AFDabs$FWH_Klet_nearest<-FayandWusH_Klet[FWH_Klet_nearestAFDabs,4]
FWH_Klet_nearestNIELSEN<-nearest(NIELSEN_GRange,FayandWusH_Klet_GRange,ignore.strand=T)
NIELSEN$FWH_Klet_nearest<-FayandWusH_Klet[FWH_Klet_nearestNIELSEN,4]
FWH_Klet_nearestVARLD<-nearest(VARLD_GRange,FayandWusH_Klet_GRange,ignore.strand=T)
VARLD$FWH_Klet_nearest<-FayandWusH_Klet[FWH_Klet_nearestVARLD,4]
FWH_Klet_nearestFLK<-nearest(FLK_GRange,FayandWusH_Klet_GRange,ignore.strand=T)
FLK$FWH_Klet_nearest<-FayandWusH_Klet[FWH_Klet_nearestFLK,4]
FWH_Klet_nearestDD<-nearest(DD_GRange,FayandWusH_Klet_GRange,ignore.strand=T)
DD$FWH_Klet_nearest<-FayandWusH_Klet[FWH_Klet_nearestDD,4]

FWH_Kowa_nearestFST<-nearest(FST_GRange,FayandWusH_Kowa_GRange,ignore.strand=T)
FST$FWH_Kowa_nearest<-FayandWusH_Kowa[FWH_Kowa_nearestFST,4]
FWH_Kowa_nearestDXY<-nearest(DXY_GRange,FayandWusH_Kowa_GRange,ignore.strand=T)
DXY$FWH_Kowa_nearest<-FayandWusH_Kowa[FWH_Kowa_nearestDXY,4]
FWH_Kowa_nearestAFDabs<-nearest(AFDabs_GRange,FayandWusH_Kowa_GRange,ignore.strand=T)
AFDabs$FWH_Kowa_nearest<-FayandWusH_Kowa[FWH_Kowa_nearestAFDabs,4]
FWH_Kowa_nearestNIELSEN<-nearest(NIELSEN_GRange,FayandWusH_Kowa_GRange,ignore.strand=T)
NIELSEN$FWH_Kowa_nearest<-FayandWusH_Kowa[FWH_Kowa_nearestNIELSEN,4]
FWH_Kowa_nearestVARLD<-nearest(VARLD_GRange,FayandWusH_Kowa_GRange,ignore.strand=T)
VARLD$FWH_Kowa_nearest<-FayandWusH_Kowa[FWH_Kowa_nearestVARLD,4]
FWH_Kowa_nearestFLK<-nearest(FLK_GRange,FayandWusH_Kowa_GRange,ignore.strand=T)
FLK$FWH_Kowa_nearest<-FayandWusH_Kowa[FWH_Kowa_nearestFLK,4]
FWH_Kowa_nearestDD<-nearest(DD_GRange,FayandWusH_Kowa_GRange,ignore.strand=T)
DD$FWH_Kowa_nearest<-FayandWusH_Kowa[FWH_Kowa_nearestDD,4]

FST$TD_diff_nearest_Klet_Kowa<-TajimasD_Klet[TD_Klet_nearestFST,4]-TajimasD_Kowa[TD_Kowa_nearestFST,4]
FST$FWH_diff_nearest_Klet_Kowa<-FayandWusH_Klet[FWH_Klet_nearestFST,4]-FayandWusH_Kowa[FWH_Kowa_nearestFST,4]
NIELSEN$TD_diff_nearest_Klet_Kowa<-TajimasD_Klet[TD_Klet_nearestNIELSEN,4]-TajimasD_Kowa[TD_Kowa_nearestNIELSEN,4]
NIELSEN$FWH_diff_nearest_Klet_Kowa<-FayandWusH_Klet[FWH_Klet_nearestNIELSEN,4]-FayandWusH_Kowa[FWH_Kowa_nearestNIELSEN,4]
DXY$TD_diff_nearest_Klet_Kowa<-TajimasD_Klet[TD_Klet_nearestDXY,4]-TajimasD_Kowa[TD_Kowa_nearestDXY,4]
DXY$FWH_diff_nearest_Klet_Kowa<-FayandWusH_Klet[FWH_Klet_nearestDXY,4]-FayandWusH_Kowa[FWH_Kowa_nearestDXY,4]
DD$TD_diff_nearest_Klet_Kowa<-TajimasD_Klet[TD_Klet_nearestDD,4]-TajimasD_Kowa[TD_Kowa_nearestDD,4]
DD$FWH_diff_nearest_Klet_Kowa<-FayandWusH_Klet[FWH_Klet_nearestDD,4]-FayandWusH_Kowa[FWH_Kowa_nearestDD,4]
VARLD$TD_diff_nearest_Klet_Kowa<-TajimasD_Klet[TD_Klet_nearestVARLD,4]-TajimasD_Kowa[TD_Kowa_nearestVARLD,4]
VARLD$FWH_diff_nearest_Klet_Kowa<-FayandWusH_Klet[FWH_Klet_nearestVARLD,4]-FayandWusH_Kowa[FWH_Kowa_nearestVARLD,4]
FLK$TD_diff_nearest_Klet_Kowa<-TajimasD_Klet[TD_Klet_nearestFLK,4]-TajimasD_Kowa[TD_Kowa_nearestFLK,4]
FLK$FWH_diff_nearest_Klet_Kowa<-FayandWusH_Klet[FWH_Klet_nearestFLK,4]-FayandWusH_Kowa[FWH_Kowa_nearestFLK,4]
AFDabs$TD_diff_nearest_Klet_Kowa<-TajimasD_Klet[TD_Klet_nearestAFDabs,4]-TajimasD_Kowa[TD_Kowa_nearestAFDabs,4]
AFDabs$FWH_diff_nearest_Klet_Kowa<-FayandWusH_Klet[FWH_Klet_nearestAFDabs,4]-FayandWusH_Kowa[FWH_Kowa_nearestAFDabs,4]



##############################
#TD and FWH standardized diff#
##############################
require(vegan)
TajimasD_Klet_stand<-decostand(TajimasD_Klet[,4],"standardize")
TajimasD_Kowa_stand<-decostand(TajimasD_Kowa[,4],"standardize")
FayandWusH_Klet_stand<-decostand(FayandWusH_Klet[,4],"standardize")
FayandWusH_Kowa_stand<-decostand(FayandWusH_Kowa[,4],"standardize")

FST<-cbind(FST,TajimasD_Klet_stand[TD_Klet_nearestFST],TajimasD_Kowa_stand[TD_Kowa_nearestFST],FayandWusH_Klet_stand[FWH_Klet_nearestFST],FayandWusH_Kowa_stand[FWH_Kowa_nearestFST])
names(FST)[43:46]<-c("TajimasD_Klet_nearest_stand","TajimasD_Kowa_nearest_stand","FayandWusH_Klet_nearest_stand","FayandWusH_Kowa_nearest_stand")
FST<-FST[,c(1:40,43:46,41:42)]
FST<-FST[,-c(36)]
NIELSEN<-cbind(NIELSEN,TajimasD_Klet_stand[TD_Klet_nearestNIELSEN],TajimasD_Kowa_stand[TD_Kowa_nearestNIELSEN],FayandWusH_Klet_stand[FWH_Klet_nearestNIELSEN],FayandWusH_Kowa_stand[FWH_Kowa_nearestNIELSEN])
names(NIELSEN)[43:46]<-c("TajimasD_Klet_nearest_stand","TajimasD_Kowa_nearest_stand","FayandWusH_Klet_nearest_stand","FayandWusH_Kowa_nearest_stand")
NIELSEN<-NIELSEN[,c(1:40,43:46,41:42)]
NIELSEN<-NIELSEN[,-c(36)]
DXY<-cbind(DXY,TajimasD_Klet_stand[TD_Klet_nearestDXY],TajimasD_Kowa_stand[TD_Kowa_nearestDXY],FayandWusH_Klet_stand[FWH_Klet_nearestDXY],FayandWusH_Kowa_stand[FWH_Kowa_nearestDXY])
names(DXY)[43:46]<-c("TajimasD_Klet_nearest_stand","TajimasD_Kowa_nearest_stand","FayandWusH_Klet_nearest_stand","FayandWusH_Kowa_nearest_stand")
DXY<-DXY[,c(1:40,43:46,41:42)]
DXY<-DXY[,-c(36)]
FLK<-cbind(FLK,TajimasD_Klet_stand[TD_Klet_nearestFLK],TajimasD_Kowa_stand[TD_Kowa_nearestFLK],FayandWusH_Klet_stand[FWH_Klet_nearestFLK],FayandWusH_Kowa_stand[FWH_Kowa_nearestFLK])
names(FLK)[43:46]<-c("TajimasD_Klet_nearest_stand","TajimasD_Kowa_nearest_stand","FayandWusH_Klet_nearest_stand","FayandWusH_Kowa_nearest_stand")
FLK<-FLK[,c(1:40,43:46,41:42)]
FLK<-FLK[,-c(36)]
DD<-cbind(DD,TajimasD_Klet_stand[TD_Klet_nearestDD],TajimasD_Kowa_stand[TD_Kowa_nearestDD],FayandWusH_Klet_stand[FWH_Klet_nearestDD],FayandWusH_Kowa_stand[FWH_Kowa_nearestDD])
names(DD)[43:46]<-c("TajimasD_Klet_nearest_stand","TajimasD_Kowa_nearest_stand","FayandWusH_Klet_nearest_stand","FayandWusH_Kowa_nearest_stand")
DD<-DD[,c(1:40,43:46,41:42)]
DD<-DD[,-c(36)]
AFDabs<-cbind(AFDabs,TajimasD_Klet_stand[TD_Klet_nearestAFDabs],TajimasD_Kowa_stand[TD_Kowa_nearestAFDabs],FayandWusH_Klet_stand[FWH_Klet_nearestAFDabs],FayandWusH_Kowa_stand[FWH_Kowa_nearestAFDabs])
names(AFDabs)[43:46]<-c("TajimasD_Klet_nearest_stand","TajimasD_Kowa_nearest_stand","FayandWusH_Klet_nearest_stand","FayandWusH_Kowa_nearest_stand")
AFDabs<-AFDabs[,c(1:40,43:46,41:42)]
AFDabs<-AFDabs[,-c(36)]
VARLD<-cbind(VARLD,TajimasD_Klet_stand[TD_Klet_nearestVARLD],TajimasD_Kowa_stand[TD_Kowa_nearestVARLD],FayandWusH_Klet_stand[FWH_Klet_nearestVARLD],FayandWusH_Kowa_stand[FWH_Kowa_nearestVARLD])
names(VARLD)[43:46]<-c("TajimasD_Klet_nearest_stand","TajimasD_Kowa_nearest_stand","FayandWusH_Klet_nearest_stand","FayandWusH_Kowa_nearest_stand")
VARLD<-VARLD[,c(1:40,43:46,41:42)]
VARLD<-VARLD[,-c(36)]

TDdiffstandFST<-TajimasD_Klet_stand[TD_Klet_nearestFST]-TajimasD_Kowa_stand[TD_Kowa_nearestFST]
FWHdiffstandFST<-FayandWusH_Klet_stand[FWH_Klet_nearestFST]-FayandWusH_Kowa_stand[FWH_Kowa_nearestFST]
TDdiffstandDXY<-TajimasD_Klet_stand[TD_Klet_nearestDXY]-TajimasD_Kowa_stand[TD_Kowa_nearestDXY]
FWHdiffstandDXY<-FayandWusH_Klet_stand[FWH_Klet_nearestDXY]-FayandWusH_Kowa_stand[FWH_Kowa_nearestDXY]
TDdiffstandAFDabs<-TajimasD_Klet_stand[TD_Klet_nearestAFDabs]-TajimasD_Kowa_stand[TD_Kowa_nearestAFDabs]
FWHdiffstandAFDabs<-FayandWusH_Klet_stand[FWH_Klet_nearestAFDabs]-FayandWusH_Kowa_stand[FWH_Kowa_nearestAFDabs]
TDdiffstandNIELSEN<-TajimasD_Klet_stand[TD_Klet_nearestNIELSEN]-TajimasD_Kowa_stand[TD_Kowa_nearestNIELSEN]
FWHdiffstandNIELSEN<-FayandWusH_Klet_stand[FWH_Klet_nearestNIELSEN]-FayandWusH_Kowa_stand[FWH_Kowa_nearestNIELSEN]
TDdiffstandVARLD<-TajimasD_Klet_stand[TD_Klet_nearestVARLD]-TajimasD_Kowa_stand[TD_Kowa_nearestVARLD]
FWHdiffstandVARLD<-FayandWusH_Klet_stand[FWH_Klet_nearestVARLD]-FayandWusH_Kowa_stand[FWH_Kowa_nearestVARLD]
TDdiffstandFLK<-TajimasD_Klet_stand[TD_Klet_nearestFLK]-TajimasD_Kowa_stand[TD_Kowa_nearestFLK]
FWHdiffstandFLK<-FayandWusH_Klet_stand[FWH_Klet_nearestFLK]-FayandWusH_Kowa_stand[FWH_Kowa_nearestFLK]
TDdiffstandDD<-TajimasD_Klet_stand[TD_Klet_nearestDD]-TajimasD_Kowa_stand[TD_Kowa_nearestDD]
FWHdiffstandDD<-FayandWusH_Klet_stand[FWH_Klet_nearestDD]-FayandWusH_Kowa_stand[FWH_Kowa_nearestDD]

FST$TD_diff_nearest_Klet_Kowa_stand<-TDdiffstandFST
FST$FWH_diff_nearest_Klet_Kowa_stand<-FWHdiffstandFST
DXY$TD_diff_nearest_Klet_Kowa_stand<-TDdiffstandDXY
DXY$FWH_diff_nearest_Klet_Kowa_stand<-FWHdiffstandDXY
AFDabs$TD_diff_nearest_Klet_Kowa_stand<-TDdiffstandAFDabs
AFDabs$FWH_diff_nearest_Klet_Kowa_stand<-FWHdiffstandAFDabs
NIELSEN$TD_diff_nearest_Klet_Kowa_stand<-TDdiffstandNIELSEN
NIELSEN$FWH_diff_nearest_Klet_Kowa_stand<-FWHdiffstandNIELSEN
VARLD$TD_diff_nearest_Klet_Kowa_stand<-TDdiffstandVARLD
VARLD$FWH_diff_nearest_Klet_Kowa_stand<-FWHdiffstandVARLD
FLK$TD_diff_nearest_Klet_Kowa_stand<-TDdiffstandFLK
FLK$FWH_diff_nearest_Klet_Kowa_stand<-FWHdiffstandFLK
DD$TD_diff_nearest_Klet_Kowa_stand<-TDdiffstandDD
DD$FWH_diff_nearest_Klet_Kowa_stand<-FWHdiffstandDD




########################################
#replace SNPeff info to window specific#
########################################
#load SNPeff data
SNPeffdf<-read.table("SNPeff_1perc_KK_arenosa.table",header=T,sep="\t")

SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]
SNPeffdfMODIFIER_GRange<-GRanges(seqnames=SNPeffdfMODIFIER$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER$POS,end=SNPeffdfMODIFIER$POS))
values(SNPeffdfMODIFIER_GRange)<-SNPeffdfMODIFIER[,1:6]

FstallMODIFIER=mergeByOverlaps(FST_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
FstallMODIFIER_df=as.data.frame(FstallMODIFIER)
FstallMODIFIER_list<-FstallMODIFIER_df[as.character(FstallMODIFIER_df$Lyr_Gene)==as.character(FstallMODIFIER_df$Gene.1),]
FST$MODIFIER<-rep("NO",nrow(FST))
FST$MODIFIER<-ifelse(!is.na(match(paste(FST$Lyr_Gene,FST$Chr,FST$Window_start,FST$Window_end),paste(FstallMODIFIER_list$Gene.1,FstallMODIFIER_list$FST_GRange.seqnames,FstallMODIFIER_list$FST_GRange.start,FstallMODIFIER_list$FST_GRange.end))),"MODIFIER",FST$MODIFIER)

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
SNPeffdfMODERATE_GRange<-GRanges(seqnames=SNPeffdfMODERATE$CHROM,ranges=IRanges(start=SNPeffdfMODERATE$POS,end=SNPeffdfMODERATE$POS))
values(SNPeffdfMODERATE_GRange)<-SNPeffdfMODERATE[,1:6]

FstallMODERATE=mergeByOverlaps(FST_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
FstallMODERATE_df=as.data.frame(FstallMODERATE)
FstallMODERATE_list<-FstallMODERATE_df[as.character(FstallMODERATE_df$Lyr_Gene)==as.character(FstallMODERATE_df$Gene.1),]
FST$MODERATE<-rep("NO",nrow(FST))
FST$MODERATE<-ifelse(!is.na(match(paste(FST$Lyr_Gene,FST$Chr,FST$Window_start,FST$Window_end),paste(FstallMODERATE_list$Gene.1,FstallMODERATE_list$FST_GRange.seqnames,FstallMODERATE_list$FST_GRange.start,FstallMODERATE_list$FST_GRange.end))),"MODERATE",FST$MODERATE)

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
SNPeffdfHIGH_GRange<-GRanges(seqnames=SNPeffdfHIGH$CHROM,ranges=IRanges(start=SNPeffdfHIGH$POS,end=SNPeffdfHIGH$POS))
values(SNPeffdfHIGH_GRange)<-SNPeffdfHIGH[,1:6]

FstallHIGH=mergeByOverlaps(FST_GRange,SNPeffdfHIGH_GRange,type=c("any"))
FstallHIGH_df=as.data.frame(FstallHIGH)
FstallHIGH_list<-FstallHIGH_df[as.character(FstallHIGH_df$Lyr_Gene)==as.character(FstallHIGH_df$Gene.1),]
FST$HIGH<-rep("NO",nrow(FST))
FST$HIGH<-ifelse(!is.na(match(paste(FST$Lyr_Gene,FST$Chr,FST$Window_start,FST$Window_end),paste(FstallHIGH_list$Gene.1,FstallHIGH_list$FST_GRange.seqnames,FstallHIGH_list$FST_GRange.start,FstallHIGH_list$FST_GRange.end))),"HIGH",FST$HIGH)

NielsenallMODIFIER=mergeByOverlaps(NIELSEN_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
NielsenallMODIFIER_df=as.data.frame(NielsenallMODIFIER)
NielsenallMODIFIER_list<-NielsenallMODIFIER_df[as.character(NielsenallMODIFIER_df$Lyr_Gene)==as.character(NielsenallMODIFIER_df$Gene.1),]
NIELSEN$MODIFIER<-rep("NO",nrow(NIELSEN))
NIELSEN$MODIFIER<-ifelse(!is.na(match(paste(NIELSEN$Lyr_Gene,NIELSEN$Chr,NIELSEN$Window_start,NIELSEN$Window_end),paste(NielsenallMODIFIER_list$Gene.1,NielsenallMODIFIER_list$NIELSEN_GRange.seqnames,NielsenallMODIFIER_list$NIELSEN_GRange.start,NielsenallMODIFIER_list$NIELSEN_GRange.end))),"MODIFIER",NIELSEN$MODIFIER)

FlkallMODIFIER=mergeByOverlaps(FLK_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
FlkallMODIFIER_df=as.data.frame(FlkallMODIFIER)
FlkallMODIFIER_list<-FlkallMODIFIER_df[as.character(FlkallMODIFIER_df$Lyr_Gene)==as.character(FlkallMODIFIER_df$Gene.1),]
FLK$MODIFIER<-rep("NO",nrow(FLK))
FLK$MODIFIER<-ifelse(!is.na(match(paste(FLK$Lyr_Gene,FLK$Chr,FLK$Window_start,FLK$Window_end),paste(FlkallMODIFIER_list$Gene.1,FlkallMODIFIER_list$FLK_GRange.seqnames,FlkallMODIFIER_list$FLK_GRange.start,FlkallMODIFIER_list$FLK_GRange.end))),"MODIFIER",FLK$MODIFIER)

DxyallMODIFIER=mergeByOverlaps(DXY_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
DxyallMODIFIER_df=as.data.frame(DxyallMODIFIER)
DxyallMODIFIER_list<-DxyallMODIFIER_df[as.character(DxyallMODIFIER_df$Lyr_Gene)==as.character(DxyallMODIFIER_df$Gene.1),]
DXY$MODIFIER<-rep("NO",nrow(DXY))
DXY$MODIFIER<-ifelse(!is.na(match(paste(DXY$Lyr_Gene,DXY$Chr,DXY$Window_start,DXY$Window_end),paste(DxyallMODIFIER_list$Gene.1,DxyallMODIFIER_list$DXY_GRange.seqnames,DxyallMODIFIER_list$DXY_GRange.start,DxyallMODIFIER_list$DXY_GRange.end))),"MODIFIER",DXY$MODIFIER)

VARLDallMODIFIER=mergeByOverlaps(VARLD_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
VARLDallMODIFIER_df=as.data.frame(VARLDallMODIFIER)
VARLDallMODIFIER_list<-VARLDallMODIFIER_df[as.character(VARLDallMODIFIER_df$Lyr_Gene)==as.character(VARLDallMODIFIER_df$Gene.1),]
VARLD$MODIFIER<-rep("NO",nrow(VARLD))
VARLD$MODIFIER<-ifelse(!is.na(match(paste(VARLD$Lyr_Gene,VARLD$Chr,VARLD$Window_start,VARLD$Window_end),paste(VARLDallMODIFIER_list$Gene.1,VARLDallMODIFIER_list$VARLD_GRange.seqnames,VARLDallMODIFIER_list$VARLD_GRange.start,VARLDallMODIFIER_list$VARLD_GRange.end))),"MODIFIER",VARLD$MODIFIER)

AFDabsallMODIFIER=mergeByOverlaps(AFDabs_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
AFDabsallMODIFIER_df=as.data.frame(AFDabsallMODIFIER)
AFDabsallMODIFIER_list<-AFDabsallMODIFIER_df[as.character(AFDabsallMODIFIER_df$Lyr_Gene)==as.character(AFDabsallMODIFIER_df$Gene.1),]
AFDabs$MODIFIER<-rep("NO",nrow(AFDabs))
AFDabs$MODIFIER<-ifelse(!is.na(match(paste(AFDabs$Lyr_Gene,AFDabs$Chr,AFDabs$Window_start,AFDabs$Window_end),paste(AFDabsallMODIFIER_list$Gene.1,AFDabsallMODIFIER_list$AFDabs_GRange.seqnames,AFDabsallMODIFIER_list$AFDabs_GRange.start,AFDabsallMODIFIER_list$AFDabs_GRange.end))),"MODIFIER",AFDabs$MODIFIER)

DDallMODIFIER=mergeByOverlaps(DD_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
DDallMODIFIER_df=as.data.frame(DDallMODIFIER)
DDallMODIFIER_list<-DDallMODIFIER_df[as.character(DDallMODIFIER_df$Lyr_Gene)==as.character(DDallMODIFIER_df$Gene.1),]
DD$MODIFIER<-rep("NO",nrow(DD))
DD$MODIFIER<-ifelse(!is.na(match(paste(DD$Lyr_Gene,DD$Chr,DD$Window_start,DD$Window_end),paste(DDallMODIFIER_list$Gene.1,DDallMODIFIER_list$DD_GRange.seqnames,DDallMODIFIER_list$DD_GRange.start,DDallMODIFIER_list$DD_GRange.end))),"MODIFIER",DD$MODIFIER)

NielsenallMODERATE=mergeByOverlaps(NIELSEN_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
NielsenallMODERATE_df=as.data.frame(NielsenallMODERATE)
NielsenallMODERATE_list<-NielsenallMODERATE_df[as.character(NielsenallMODERATE_df$Lyr_Gene)==as.character(NielsenallMODERATE_df$Gene.1),]
NIELSEN$MODERATE<-rep("NO",nrow(NIELSEN))
NIELSEN$MODERATE<-ifelse(match(paste(NIELSEN$Lyr_Gene,NIELSEN$Chr,NIELSEN$Window_start,NIELSEN$Window_end),paste(NielsenallMODERATE_list$Gene.1,NielsenallMODERATE_list$NIELSEN_GRange.seqnames,NielsenallMODERATE_list$NIELSEN_GRange.start,NielsenallMODERATE_list$NIELSEN_GRange.end)),"MODERATE",NIELSEN$MODERATE)

FlkallMODERATE=mergeByOverlaps(FLK_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
FlkallMODERATE_df=as.data.frame(FlkallMODERATE)
FlkallMODERATE_list<-FlkallMODERATE_df[as.character(FlkallMODERATE_df$Lyr_Gene)==as.character(FlkallMODERATE_df$Gene.1),]
FLK$MODERATE<-rep("NO",nrow(FLK))
FLK$MODERATE<-ifelse(!is.na(match(paste(FLK$Lyr_Gene,FLK$Chr,FLK$Window_start,FLK$Window_end),paste(FlkallMODERATE_list$Gene.1,FlkallMODERATE_list$FLK_GRange.seqnames,FlkallMODERATE_list$FLK_GRange.start,FlkallMODERATE_list$FLK_GRange.end))),"MODERATE",FLK$MODERATE)

DxyallMODERATE=mergeByOverlaps(DXY_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
DxyallMODERATE_df=as.data.frame(DxyallMODERATE)
DxyallMODERATE_list<-DxyallMODERATE_df[as.character(DxyallMODERATE_df$Lyr_Gene)==as.character(DxyallMODERATE_df$Gene.1),]
DXY$MODERATE<-rep("NO",nrow(DXY))
DXY$MODERATE<-ifelse(!is.na(match(paste(DXY$Lyr_Gene,DXY$Chr,DXY$Window_start,DXY$Window_end),paste(DxyallMODERATE_list$Gene.1,DxyallMODERATE_list$DXY_GRange.seqnames,DxyallMODERATE_list$DXY_GRange.start,DxyallMODERATE_list$DXY_GRange.end))),"MODERATE",DXY$MODERATE)

VARLDallMODERATE=mergeByOverlaps(VARLD_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
VARLDallMODERATE_df=as.data.frame(VARLDallMODERATE)
VARLDallMODERATE_list<-VARLDallMODERATE_df[as.character(VARLDallMODERATE_df$Lyr_Gene)==as.character(VARLDallMODERATE_df$Gene.1),]
VARLD$MODERATE<-rep("NO",nrow(VARLD))
VARLD$MODERATE<-ifelse(!is.na(match(paste(VARLD$Lyr_Gene,VARLD$Chr,VARLD$Window_start,VARLD$Window_end),paste(VARLDallMODERATE_list$Gene.1,VARLDallMODERATE_list$VARLD_GRange.seqnames,VARLDallMODERATE_list$VARLD_GRange.start,VARLDallMODERATE_list$VARLD_GRange.end))),"MODERATE",VARLD$MODERATE)

AFDabsallMODERATE=mergeByOverlaps(AFDabs_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
AFDabsallMODERATE_df=as.data.frame(AFDabsallMODERATE)
AFDabsallMODERATE_list<-AFDabsallMODERATE_df[as.character(AFDabsallMODERATE_df$Lyr_Gene)==as.character(AFDabsallMODERATE_df$Gene.1),]
AFDabs$MODERATE<-rep("NO",nrow(AFDabs))
AFDabs$MODERATE<-ifelse(!is.na(match(paste(AFDabs$Lyr_Gene,AFDabs$Chr,AFDabs$Window_start,AFDabs$Window_end),paste(AFDabsallMODERATE_list$Gene.1,AFDabsallMODERATE_list$AFDabs_GRange.seqnames,AFDabsallMODERATE_list$AFDabs_GRange.start,AFDabsallMODERATE_list$AFDabs_GRange.end))),"MODERATE",AFDabs$MODERATE)

DDallMODERATE=mergeByOverlaps(DD_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
DDallMODERATE_df=as.data.frame(DDallMODERATE)
DDallMODERATE_list<-DDallMODERATE_df[as.character(DDallMODERATE_df$Lyr_Gene)==as.character(DDallMODERATE_df$Gene.1),]
DD$MODERATE<-rep("NO",nrow(DD))
DD$MODERATE<-ifelse(!is.na(match(paste(DD$Lyr_Gene,DD$Chr,DD$Window_start,DD$Window_end),paste(DDallMODERATE_list$Gene.1,DDallMODERATE_list$DD_GRange.seqnames,DDallMODERATE_list$DD_GRange.start,DDallMODERATE_list$DD_GRange.end))),"MODERATE",DD$MODERATE)

NielsenallHIGH=mergeByOverlaps(NIELSEN_GRange,SNPeffdfHIGH_GRange,type=c("any"))
NielsenallHIGH_df=as.data.frame(NielsenallHIGH)
NielsenallHIGH_list<-NielsenallHIGH_df[as.character(NielsenallHIGH_df$Lyr_Gene)==as.character(NielsenallHIGH_df$Gene.1),]
NIELSEN$HIGH<-rep("NO",nrow(NIELSEN))
NIELSEN$HIGH<-ifelse(!is.na(match(paste(NIELSEN$Lyr_Gene,NIELSEN$Chr,NIELSEN$Window_start,NIELSEN$Window_end),paste(NielsenallHIGH_list$Gene.1,NielsenallHIGH_list$NIELSEN_GRange.seqnames,NielsenallHIGH_list$NIELSEN_GRange.start,NielsenallHIGH_list$NIELSEN_GRange.end))),"HIGH",NIELSEN$HIGH)

FlkallHIGH=mergeByOverlaps(FLK_GRange,SNPeffdfHIGH_GRange,type=c("any"))
FlkallHIGH_df=as.data.frame(FlkallHIGH)
FlkallHIGH_list<-FlkallHIGH_df[as.character(FlkallHIGH_df$Lyr_Gene)==as.character(FlkallHIGH_df$Gene.1),]
FLK$HIGH<-rep("NO",nrow(FLK))
FLK$HIGH<-ifelse(!is.na(match(paste(FLK$Lyr_Gene,FLK$Chr,FLK$Window_start,FLK$Window_end),paste(FlkallHIGH_list$Gene.1,FlkallHIGH_list$FLK_GRange.seqnames,FlkallHIGH_list$FLK_GRange.start,FlkallHIGH_list$FLK_GRange.end))),"HIGH",FLK$HIGH)

DxyallHIGH=mergeByOverlaps(DXY_GRange,SNPeffdfHIGH_GRange,type=c("any"))
DxyallHIGH_df=as.data.frame(DxyallHIGH)
DxyallHIGH_list<-DxyallHIGH_df[as.character(DxyallHIGH_df$Lyr_Gene)==as.character(DxyallHIGH_df$Gene.1),]
DXY$HIGH<-rep("NO",nrow(DXY))
DXY$HIGH<-ifelse(!is.na(match(paste(DXY$Lyr_Gene,DXY$Chr,DXY$Window_start,DXY$Window_end),paste(DxyallHIGH_list$Gene.1,DxyallHIGH_list$DXY_GRange.seqnames,DxyallHIGH_list$DXY_GRange.start,DxyallHIGH_list$DXY_GRange.end))),"HIGH",DXY$HIGH)

VARLDallHIGH=mergeByOverlaps(VARLD_GRange,SNPeffdfHIGH_GRange,type=c("any"))
VARLDallHIGH_df=as.data.frame(VARLDallHIGH)
VARLDallHIGH_list<-VARLDallHIGH_df[as.character(VARLDallHIGH_df$Lyr_Gene)==as.character(VARLDallHIGH_df$Gene.1),]
VARLD$HIGH<-rep("NO",nrow(VARLD))
VARLD$HIGH<-ifelse(!is.na(match(paste(VARLD$Lyr_Gene,VARLD$Chr,VARLD$Window_start,VARLD$Window_end),paste(VARLDallHIGH_list$Gene.1,VARLDallHIGH_list$VARLD_GRange.seqnames,VARLDallHIGH_list$VARLD_GRange.start,VARLDallHIGH_list$VARLD_GRange.end))),"HIGH",VARLD$HIGH)

AFDabsallHIGH=mergeByOverlaps(AFDabs_GRange,SNPeffdfHIGH_GRange,type=c("any"))
AFDabsallHIGH_df=as.data.frame(AFDabsallHIGH)
AFDabsallHIGH_list<-AFDabsallHIGH_df[as.character(AFDabsallHIGH_df$Lyr_Gene)==as.character(AFDabsallHIGH_df$Gene.1),]
AFDabs$HIGH<-rep("NO",nrow(AFDabs))
AFDabs$HIGH<-ifelse(!is.na(match(paste(AFDabs$Lyr_Gene,AFDabs$Chr,AFDabs$Window_start,AFDabs$Window_end),paste(AFDabsallHIGH_list$Gene.1,AFDabsallHIGH_list$AFDabs_GRange.seqnames,AFDabsallHIGH_list$AFDabs_GRange.start,AFDabsallHIGH_list$AFDabs_GRange.end))),"HIGH",AFDabs$HIGH)

DDallHIGH=mergeByOverlaps(DD_GRange,SNPeffdfHIGH_GRange,type=c("any"))
DDallHIGH_df=as.data.frame(DDallHIGH)
DDallHIGH_list<-DDallHIGH_df[as.character(DDallHIGH_df$Lyr_Gene)==as.character(DDallHIGH_df$Gene.1),]
DD$HIGH<-rep("NO",nrow(DD))
DD$HIGH<-ifelse(!is.na(match(paste(DD$Lyr_Gene,DD$Chr,DD$Window_start,DD$Window_end),paste(DDallHIGH_list$Gene.1,DDallHIGH_list$DD_GRange.seqnames,DDallHIGH_list$DD_GRange.start,DDallHIGH_list$DD_GRange.end))),"HIGH",DD$HIGH)


#############################################################################
#add SNPeff info with info for whole gene#
#############################################################################

n1<-36 #number of Klet individuals * ploidy
n2<-32 #number of Kowa individuals * ploidy
AFdata<-read.table("LyrataKletKowaGS_woGQfil_arenosa.csv",header=TRUE,sep="\t")
test<-as.character(AFdata$Chrom)
AFdata[,1]<-test
AFdata[,3]<-AFdata[,3]/n1
AFdata[,6]<-AFdata[,6]/n2
AFdata2<-data.frame(AFdata[,1:3],AFdata[,6])
names(AFdata2)<-c("CHROM","POS","AF_Klet","AF_Kowa")

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
SNPeffdfHIGH2<-merge(SNPeffdfHIGH,AFdata2)

FSTHighKlet<-vector(mode="character",length=nrow(FST))
for (i in 1:nrow(FST)) {
	if(FST$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		FSTHighKlet[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Klet[match(FST$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		FSTHighKlet[i]<-"NO"
		}
	}
FSTHighKowa<-vector(mode="character",length=nrow(FST))
for (i in 1:nrow(FST)) {
	if(FST$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		FSTHighKowa[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kowa[match(FST$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		FSTHighKowa[i]<-"NO"
		}
	}

test<-FST
FST<-data.frame(append(test,list(AF_HIGH_Klet_gene=FSTHighKlet),after=match("HIGH",names(test))))
FST<-data.frame(append(FST,list(AF_HIGH_Kowa_gene=FSTHighKowa),after=match("AF_HIGH_Klet_gene",names(FST))))

NIELSENHighKlet<-vector(mode="character",length=nrow(NIELSEN))
for (i in 1:nrow(NIELSEN)) {
	if(NIELSEN$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		NIELSENHighKlet[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Klet[match(NIELSEN$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		NIELSENHighKlet[i]<-"NO"
		}
	}
NIELSENHighKowa<-vector(mode="character",length=nrow(NIELSEN))
for (i in 1:nrow(NIELSEN)) {
	if(NIELSEN$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		NIELSENHighKowa[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kowa[match(NIELSEN$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		NIELSENHighKowa[i]<-"NO"
		}
	}

test<-NIELSEN
NIELSEN<-data.frame(append(test,list(AF_HIGH_Klet_gene=NIELSENHighKlet),after=match("HIGH",names(test))))
NIELSEN<-data.frame(append(NIELSEN,list(AF_HIGH_Kowa_gene=NIELSENHighKowa),after=match("AF_HIGH_Klet_gene",names(NIELSEN))))

DXYHighKlet<-vector(mode="character",length=nrow(DXY))
for (i in 1:nrow(DXY)) {
	if(DXY$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		DXYHighKlet[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Klet[match(DXY$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		DXYHighKlet[i]<-"NO"
		}
	}
DXYHighKowa<-vector(mode="character",length=nrow(DXY))
for (i in 1:nrow(DXY)) {
	if(DXY$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		DXYHighKowa[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kowa[match(DXY$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		DXYHighKowa[i]<-"NO"
		}
	}

test<-DXY
DXY<-data.frame(append(test,list(AF_HIGH_Klet_gene=DXYHighKlet),after=match("HIGH",names(test))))
DXY<-data.frame(append(DXY,list(AF_HIGH_Kowa_gene=DXYHighKowa),after=match("AF_HIGH_Klet_gene",names(DXY))))

AFDabsHighKlet<-vector(mode="character",length=nrow(AFDabs))
for (i in 1:nrow(AFDabs)) {
	if(AFDabs$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		AFDabsHighKlet[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Klet[match(AFDabs$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		AFDabsHighKlet[i]<-"NO"
		}
	}
AFDabsHighKowa<-vector(mode="character",length=nrow(AFDabs))
for (i in 1:nrow(AFDabs)) {
	if(AFDabs$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		AFDabsHighKowa[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kowa[match(AFDabs$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		AFDabsHighKowa[i]<-"NO"
		}
	}

test<-AFDabs
AFDabs<-data.frame(append(test,list(AF_HIGH_Klet_gene=AFDabsHighKlet),after=match("HIGH",names(test))))
AFDabs<-data.frame(append(AFDabs,list(AF_HIGH_Kowa_gene=AFDabsHighKowa),after=match("AF_HIGH_Klet_gene",names(AFDabs))))

FLKHighKlet<-vector(mode="character",length=nrow(FLK))
for (i in 1:nrow(FLK)) {
	if(FLK$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		FLKHighKlet[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Klet[match(FLK$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		FLKHighKlet[i]<-"NO"
		}
	}
FLKHighKowa<-vector(mode="character",length=nrow(FLK))
for (i in 1:nrow(FLK)) {
	if(FLK$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		FLKHighKowa[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kowa[match(FLK$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		FLKHighKowa[i]<-"NO"
		}
	}

test<-FLK
FLK<-data.frame(append(test,list(AF_HIGH_Klet_gene=FLKHighKlet),after=match("HIGH",names(test))))
FLK<-data.frame(append(FLK,list(AF_HIGH_Kowa_gene=FLKHighKowa),after=match("AF_HIGH_Klet_gene",names(FLK))))

DDHighKlet<-vector(mode="character",length=nrow(DD))
for (i in 1:nrow(DD)) {
	if(DD$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		DDHighKlet[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Klet[match(DD$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		DDHighKlet[i]<-"NO"
		}
	}
DDHighKowa<-vector(mode="character",length=nrow(DD))
for (i in 1:nrow(DD)) {
	if(DD$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		DDHighKowa[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kowa[match(DD$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		DDHighKowa[i]<-"NO"
		}
	}

test<-DD
DD<-data.frame(append(test,list(AF_HIGH_Klet_gene=DDHighKlet),after=match("HIGH",names(test))))
DD<-data.frame(append(DD,list(AF_HIGH_Kowa_gene=DDHighKowa),after=match("AF_HIGH_Klet_gene",names(DD))))

VARLDHighKlet<-vector(mode="character",length=nrow(VARLD))
for (i in 1:nrow(VARLD)) {
	if(VARLD$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		VARLDHighKlet[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Klet[match(VARLD$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		VARLDHighKlet[i]<-"NO"
		}
	}
VARLDHighKowa<-vector(mode="character",length=nrow(VARLD))
for (i in 1:nrow(VARLD)) {
	if(VARLD$Lyr_Gene[i]%in%SNPeffdfHIGH2$Gene) {
		VARLDHighKowa[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kowa[match(VARLD$Lyr_Gene[i],SNPeffdfHIGH2$Gene)]))
		} else {
		VARLDHighKowa[i]<-"NO"
		}
	}

test<-VARLD
VARLD<-data.frame(append(test,list(AF_HIGH_Klet_gene=VARLDHighKlet),after=match("HIGH",names(test))))
VARLD<-data.frame(append(VARLD,list(AF_HIGH_Kowa_gene=VARLDHighKowa),after=match("AF_HIGH_Klet_gene",names(VARLD))))


#############################################
#add Rank info#
#############################################
FST$Rank_Fst<-rank(-as.numeric(as.character(FST$Fst)),ties.method="min")
NIELSEN$Rank_Nielsen<-rank(-as.numeric(as.character(NIELSEN$Nielsen)),ties.method="min")
DXY$Rank_Dxy<-rank(-as.numeric(as.character(DXY$Dxy)),ties.method="min")
DD$Rank_DD<-rank(as.numeric(as.character(DD$DD)),ties.method="min")
AFDabs$Rank_AFDabs<-rank(-as.numeric(as.character(AFDabs$AFDabs)),ties.method="min")
FLK$Rank_Flk<-rank(-as.numeric(as.character(FLK$Flk)),ties.method="min")
VARLD$Rank_Varld<-rank(-as.numeric(as.character(VARLD$VarLD)),ties.method="min")

######################################################
#add Sweedinfo#
######################################################
dataSweed_Klet1<-read.table("Klet_Sweed_scaffold1.table",sep="\t",header=T)
dataSweed_Klet1$Scaffold<-rep("Scaffold_1",nrow(dataSweed_Klet1))
dataSweed_Klet2<-read.table("Klet_Sweed_scaffold2.table",sep="\t",header=T)
dataSweed_Klet2$Scaffold<-rep("Scaffold_2",nrow(dataSweed_Klet2))
dataSweed_Klet3<-read.table("Klet_Sweed_scaffold3.table",sep="\t",header=T)
dataSweed_Klet3$Scaffold<-rep("Scaffold_3",nrow(dataSweed_Klet3))
dataSweed_Klet4<-read.table("Klet_Sweed_scaffold4.table",sep="\t",header=T)
dataSweed_Klet4$Scaffold<-rep("Scaffold_4",nrow(dataSweed_Klet4))
dataSweed_Klet5<-read.table("Klet_Sweed_scaffold5.table",sep="\t",header=T)
dataSweed_Klet5$Scaffold<-rep("Scaffold_5",nrow(dataSweed_Klet5))
dataSweed_Klet6<-read.table("Klet_Sweed_scaffold6.table",sep="\t",header=T)
dataSweed_Klet6$Scaffold<-rep("Scaffold_6",nrow(dataSweed_Klet6))
dataSweed_Klet7<-read.table("Klet_Sweed_scaffold7.table",sep="\t",header=T)
dataSweed_Klet7$Scaffold<-rep("Scaffold_7",nrow(dataSweed_Klet7))
dataSweed_Klet8<-read.table("Klet_Sweed_scaffold8.table",sep="\t",header=T)
dataSweed_Klet8$Scaffold<-rep("Scaffold_8",nrow(dataSweed_Klet8))

dataSweed_Klet<-rbind(dataSweed_Klet1,dataSweed_Klet2,dataSweed_Klet3,dataSweed_Klet4,dataSweed_Klet5,dataSweed_Klet6,dataSweed_Klet7,dataSweed_Klet8)
dataSweed_Klet<-dataSweed_Klet[,c(4,1:3)]
dataSweed_Klet[,1]<-as.factor(dataSweed_Klet[,1])

Sweed_Klet_GR<-GRanges(seqnames=tolower(dataSweed_Klet$Scaffold),ranges=IRanges(start=dataSweed_Klet$Position,end=dataSweed_Klet$Position))
values(Sweed_Klet_GR)<-dataSweed_Klet[,3:4]

Sweed_Klet_nearestFST<-nearest(FST_GRange,Sweed_Klet_GR,ignore.strand=T)
FST$Sweed_Klet_nearest_Likelihood<-dataSweed_Klet[Sweed_Klet_nearestFST,3]
FST$Sweed_Klet_nearest_alpha<-dataSweed_Klet[Sweed_Klet_nearestFST,4]
Sweed_Klet_nearestDD<-nearest(DD_GRange,Sweed_Klet_GR,ignore.strand=T)
DD$Sweed_Klet_nearest_Likelihood<-dataSweed_Klet[Sweed_Klet_nearestDD,3]
DD$Sweed_Klet_nearest_alpha<-dataSweed_Klet[Sweed_Klet_nearestDD,4]
Sweed_Klet_nearestDXY<-nearest(DXY_GRange,Sweed_Klet_GR,ignore.strand=T)
DXY$Sweed_Klet_nearest_Likelihood<-dataSweed_Klet[Sweed_Klet_nearestDXY,3]
DXY$Sweed_Klet_nearest_alpha<-dataSweed_Klet[Sweed_Klet_nearestDXY,4]
Sweed_Klet_nearestNIELSEN<-nearest(NIELSEN_GRange,Sweed_Klet_GR,ignore.strand=T)
NIELSEN$Sweed_Klet_nearest_Likelihood<-dataSweed_Klet[Sweed_Klet_nearestNIELSEN,3]
NIELSEN$Sweed_Klet_nearest_alpha<-dataSweed_Klet[Sweed_Klet_nearestNIELSEN,4]
Sweed_Klet_nearestFLK<-nearest(FLK_GRange,Sweed_Klet_GR,ignore.strand=T)
FLK$Sweed_Klet_nearest_Likelihood<-dataSweed_Klet[Sweed_Klet_nearestFLK,3]
FLK$Sweed_Klet_nearest_alpha<-dataSweed_Klet[Sweed_Klet_nearestFLK,4]
Sweed_Klet_nearestVARLD<-nearest(VARLD_GRange,Sweed_Klet_GR,ignore.strand=T)
VARLD$Sweed_Klet_nearest_Likelihood<-dataSweed_Klet[Sweed_Klet_nearestVARLD,3]
VARLD$Sweed_Klet_nearest_alpha<-dataSweed_Klet[Sweed_Klet_nearestVARLD,4]
Sweed_Klet_nearestAFDabs<-nearest(AFDabs_GRange,Sweed_Klet_GR,ignore.strand=T)
AFDabs$Sweed_Klet_nearest_Likelihood<-dataSweed_Klet[Sweed_Klet_nearestAFDabs,3]
AFDabs$Sweed_Klet_nearest_alpha<-dataSweed_Klet[Sweed_Klet_nearestAFDabs,4]

dataSweed_Kowa1<-read.table("Kowa_Sweed_scaffold1.table",sep="\t",header=T)
dataSweed_Kowa1$Scaffold<-rep("Scaffold_1",nrow(dataSweed_Kowa1))
dataSweed_Kowa2<-read.table("Kowa_Sweed_scaffold2.table",sep="\t",header=T)
dataSweed_Kowa2$Scaffold<-rep("Scaffold_2",nrow(dataSweed_Kowa2))
dataSweed_Kowa3<-read.table("Kowa_Sweed_scaffold3.table",sep="\t",header=T)
dataSweed_Kowa3$Scaffold<-rep("Scaffold_3",nrow(dataSweed_Kowa3))
dataSweed_Kowa4<-read.table("Kowa_Sweed_scaffold4.table",sep="\t",header=T)
dataSweed_Kowa4$Scaffold<-rep("Scaffold_4",nrow(dataSweed_Kowa4))
dataSweed_Kowa5<-read.table("Kowa_Sweed_scaffold5.table",sep="\t",header=T)
dataSweed_Kowa5$Scaffold<-rep("Scaffold_5",nrow(dataSweed_Kowa5))
dataSweed_Kowa6<-read.table("Kowa_Sweed_scaffold6.table",sep="\t",header=T)
dataSweed_Kowa6$Scaffold<-rep("Scaffold_6",nrow(dataSweed_Kowa6))
dataSweed_Kowa7<-read.table("Kowa_Sweed_scaffold7.table",sep="\t",header=T)
dataSweed_Kowa7$Scaffold<-rep("Scaffold_7",nrow(dataSweed_Kowa7))
dataSweed_Kowa8<-read.table("Kowa_Sweed_scaffold8.table",sep="\t",header=T)
dataSweed_Kowa8$Scaffold<-rep("Scaffold_8",nrow(dataSweed_Kowa8))

dataSweed_Kowa<-rbind(dataSweed_Kowa1,dataSweed_Kowa2,dataSweed_Kowa3,dataSweed_Kowa4,dataSweed_Kowa5,dataSweed_Kowa6,dataSweed_Kowa7,dataSweed_Kowa8)
dataSweed_Kowa<-dataSweed_Kowa[,c(4,1:3)]
dataSweed_Kowa[,1]<-as.factor(dataSweed_Kowa[,1])

Sweed_Kowa_GR<-GRanges(seqnames=tolower(dataSweed_Kowa$Scaffold),ranges=IRanges(start=dataSweed_Kowa$Position,end=dataSweed_Kowa$Position))
values(Sweed_Kowa_GR)<-dataSweed_Kowa[,3:4]

Sweed_Kowa_nearestFST<-nearest(FST_GRange,Sweed_Kowa_GR,ignore.strand=T)
FST$Sweed_Kowa_nearest_Likelihood<-dataSweed_Kowa[Sweed_Kowa_nearestFST,3]
FST$Sweed_Kowa_nearest_alpha<-dataSweed_Kowa[Sweed_Kowa_nearestFST,4]
Sweed_Kowa_nearestDD<-nearest(DD_GRange,Sweed_Kowa_GR,ignore.strand=T)
DD$Sweed_Kowa_nearest_Likelihood<-dataSweed_Kowa[Sweed_Kowa_nearestDD,3]
DD$Sweed_Kowa_nearest_alpha<-dataSweed_Kowa[Sweed_Kowa_nearestDD,4]
Sweed_Kowa_nearestDXY<-nearest(DXY_GRange,Sweed_Kowa_GR,ignore.strand=T)
DXY$Sweed_Kowa_nearest_Likelihood<-dataSweed_Kowa[Sweed_Kowa_nearestDXY,3]
DXY$Sweed_Kowa_nearest_alpha<-dataSweed_Kowa[Sweed_Kowa_nearestDXY,4]
Sweed_Kowa_nearestNIELSEN<-nearest(NIELSEN_GRange,Sweed_Kowa_GR,ignore.strand=T)
NIELSEN$Sweed_Kowa_nearest_Likelihood<-dataSweed_Kowa[Sweed_Kowa_nearestNIELSEN,3]
NIELSEN$Sweed_Kowa_nearest_alpha<-dataSweed_Kowa[Sweed_Kowa_nearestNIELSEN,4]
Sweed_Kowa_nearestFLK<-nearest(FLK_GRange,Sweed_Kowa_GR,ignore.strand=T)
FLK$Sweed_Kowa_nearest_Likelihood<-dataSweed_Kowa[Sweed_Kowa_nearestFLK,3]
FLK$Sweed_Kowa_nearest_alpha<-dataSweed_Kowa[Sweed_Kowa_nearestFLK,4]
Sweed_Kowa_nearestVARLD<-nearest(VARLD_GRange,Sweed_Kowa_GR,ignore.strand=T)
VARLD$Sweed_Kowa_nearest_Likelihood<-dataSweed_Kowa[Sweed_Kowa_nearestVARLD,3]
VARLD$Sweed_Kowa_nearest_alpha<-dataSweed_Kowa[Sweed_Kowa_nearestVARLD,4]
Sweed_Kowa_nearestAFDabs<-nearest(AFDabs_GRange,Sweed_Kowa_GR,ignore.strand=T)
AFDabs$Sweed_Kowa_nearest_Likelihood<-dataSweed_Kowa[Sweed_Kowa_nearestAFDabs,3]
AFDabs$Sweed_Kowa_nearest_alpha<-dataSweed_Kowa[Sweed_Kowa_nearestAFDabs,4]


#################################################
#add info in how many metrics among 0.1% outlier#
#################################################

countFST<-vector(mode="numeric",length=nrow(FST))
for (i in 1:nrow(FST))
	{countFST[i]=0
	if (as.numeric(as.character(FST$Nielsen[i]))>=243.271265489257)
		{countFST[i]=countFST[i]+1
		}
	if (as.numeric(as.character(FST$AFDabs[i]))>=0.288629444444442)
		{countFST[i]=countFST[i]+1
		}
	if (as.numeric(as.character(FST$Dxy[i]))>=0.419311666666666)
		{countFST[i]=countFST[i]+1
		}
	if (as.numeric(as.character(FST$VarLD[i]))>=34.4780176)
		{countFST[i]=countFST[i]+1
		}
	if (as.numeric(as.character(FST$Flk[i]))>=3.29076116503269)
		{countFST[i]=countFST[i]+1
		}
	if(as.numeric(as.character(FST$DD[i]))<=-0.147369851025826)
		{countFST[i]=countFST[i]+1
		}
	}
FST$Number_of_metrics_0.1percent<-countFST


countNIELSEN<-vector(mode="numeric",length=nrow(NIELSEN))
for (i in 1:nrow(NIELSEN))
	{countNIELSEN[i]=0
	if (as.numeric(as.character(NIELSEN$Fst[i]))>=0.159584258881707)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	if (as.numeric(as.character(NIELSEN$AFDabs[i]))>=0.288629444444442)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	if (as.numeric(as.character(NIELSEN$Dxy[i]))>=0.419311666666666)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	if (as.numeric(as.character(NIELSEN$VarLD[i]))>=34.4780176)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	if (as.numeric(as.character(NIELSEN$Flk[i]))>=3.29076116503269)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	if(as.numeric(as.character(NIELSEN$DD[i]))<=-0.147369851025826)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	}
NIELSEN$Number_of_metrics_0.1percent<-countNIELSEN

countDXY<-vector(mode="numeric",length=nrow(DXY))
for (i in 1:nrow(DXY))
	{countDXY[i]=0
	if (as.numeric(as.character(DXY$Fst[i]))>=0.159584258881707)
		{countDXY[i]=countDXY[i]+1
		}
	if (as.numeric(as.character(DXY$AFDabs[i]))>=0.288629444444442)
		{countDXY[i]=countDXY[i]+1
		}
	if (as.numeric(as.character(DXY$Nielsen[i]))>=243.271265489257)
		{countDXY[i]=countDXY[i]+1
		}
	if (as.numeric(as.character(DXY$VarLD[i]))>=34.4780176)
		{countDXY[i]=countDXY[i]+1
		}
	if (as.numeric(as.character(DXY$Flk[i]))>=3.29076116503269)
		{countDXY[i]=countDXY[i]+1
		}
	if(as.numeric(as.character(DXY$DD[i]))<=-0.147369851025826)
		{countDXY[i]=countDXY[i]+1
		}
	}
DXY$Number_of_metrics_0.1percent<-countDXY

countVARLD<-vector(mode="numeric",length=nrow(VARLD))
for (i in 1:nrow(VARLD))
	{countVARLD[i]=0
	if (as.numeric(as.character(VARLD$Fst[i]))>=0.159584258881707)
		{countVARLD[i]=countVARLD[i]+1
		}
	if (as.numeric(as.character(VARLD$AFDabs[i]))>=0.288629444444442)
		{countVARLD[i]=countVARLD[i]+1
		}
	if (as.numeric(as.character(VARLD$Nielsen[i]))>=243.271265489257)
		{countVARLD[i]=countVARLD[i]+1
		}
	if (as.numeric(as.character(VARLD$Dxy[i]))>=0.419311666666666)
		{countVARLD[i]=countVARLD[i]+1
		}
	if (as.numeric(as.character(VARLD$Flk[i]))>=3.29076116503269)
		{countVARLD[i]=countVARLD[i]+1
		}
	if(as.numeric(as.character(VARLD$DD[i]))<=-0.147369851025826)
		{countVARLD[i]=countVARLD[i]+1
		}
	}
VARLD$Number_of_metrics_0.1percent<-countVARLD

countFLK<-vector(mode="numeric",length=nrow(FLK))
for (i in 1:nrow(FLK))
	{countFLK[i]=0
	if (as.numeric(as.character(FLK$Fst[i]))>=0.159584258881707)
		{countFLK[i]=countFLK[i]+1
		}
	if (as.numeric(as.character(FLK$AFDabs[i]))>=0.288629444444442)
		{countFLK[i]=countFLK[i]+1
		}
	if (as.numeric(as.character(FLK$Nielsen[i]))>=243.271265489257)
		{countFLK[i]=countFLK[i]+1
		}
	if (as.numeric(as.character(FLK$Dxy[i]))>=0.419311666666666)
		{countFLK[i]=countFLK[i]+1
		}
	if (as.numeric(as.character(FLK$VarLD[i]))>=34.4780176)
		{countFLK[i]=countFLK[i]+1
		}
	if(as.numeric(as.character(FLK$DD[i]))<=-0.147369851025826)
		{countFLK[i]=countFLK[i]+1
		}
	}
FLK$Number_of_metrics_0.1percent<-countFLK

countDD<-vector(mode="numeric",length=nrow(DD))
for (i in 1:nrow(DD))
	{countDD[i]=0
	if (as.numeric(as.character(DD$Fst[i]))>=0.159584258881707)
		{countDD[i]=countDD[i]+1
		}
	if (as.numeric(as.character(DD$AFDabs[i]))>=0.288629444444442)
		{countDD[i]=countDD[i]+1
		}
	if (as.numeric(as.character(DD$Nielsen[i]))>=243.271265489257)
		{countDD[i]=countDD[i]+1
		}
	if (as.numeric(as.character(DD$Dxy[i]))>=0.419311666666666)
		{countDD[i]=countDD[i]+1
		}
	if (as.numeric(as.character(DD$VarLD[i]))>=34.4780176)
		{countDD[i]=countDD[i]+1
		}
	if(as.numeric(as.character(DD$Flk[i]))>=3.29076116503269)
		{countDD[i]=countDD[i]+1
		}
	}
DD$Number_of_metrics_0.1percent<-countDD


countAFDabs<-vector(mode="numeric",length=nrow(AFDabs))
for (i in 1:nrow(AFDabs))
	{countAFDabs[i]=0
	if (as.numeric(as.character(AFDabs$Fst[i]))>=0.159584258881707)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	if (as.numeric(as.character(AFDabs$Flk[i]))>=3.29076116503269)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	if (as.numeric(as.character(AFDabs$Nielsen[i]))>=243.271265489257)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	if (as.numeric(as.character(AFDabs$Dxy[i]))>=0.419311666666666)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	if (as.numeric(as.character(AFDabs$VarLD[i]))>=34.4780176)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	if(as.numeric(as.character(AFDabs$DD[i]))<=-0.147369851025826)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	}
AFDabs$Number_of_metrics_0.1percent<-countAFDabs

#########################################
#add metal homeostasis gene list info#
#########################################


Methomlist<-read.xlsx2("metal_homeostasis_Ute_16_01_2017late.xlsx",1,header=T)

FST$Gene<-toupper(FST$Gene)
Methomlist$AGI_Number<-toupper(Methomlist$AGI_Number)
FSTmetlist<-merge(FST,Methomlist,all.x=T,by.x="Gene",by.y="AGI_Number")

NIELSEN$Gene<-toupper(NIELSEN$Gene)
Methomlist$AGI_Number<-toupper(Methomlist$AGI_Number)
NIELSENmetlist<-merge(NIELSEN,Methomlist,all.x=T,by.x="Gene",by.y="AGI_Number")

DXY$Gene<-toupper(DXY$Gene)
Methomlist$AGI_Number<-toupper(Methomlist$AGI_Number)
DXYmetlist<-merge(DXY,Methomlist,all.x=T,by.x="Gene",by.y="AGI_Number")

FLK$Gene<-toupper(FLK$Gene)
Methomlist$AGI_Number<-toupper(Methomlist$AGI_Number)
FLKmetlist<-merge(FLK,Methomlist,all.x=T,by.x="Gene",by.y="AGI_Number")

VARLD$Gene<-toupper(VARLD$Gene)
Methomlist$AGI_Number<-toupper(Methomlist$AGI_Number)
VARLDmetlist<-merge(VARLD,Methomlist,all.x=T,by.x="Gene",by.y="AGI_Number")

AFDabs$Gene<-toupper(AFDabs$Gene)
Methomlist$AGI_Number<-toupper(Methomlist$AGI_Number)
AFDabsmetlist<-merge(AFDabs,Methomlist,all.x=T,by.x="Gene",by.y="AGI_Number")

DD$Gene<-toupper(DD$Gene)
Methomlist$AGI_Number<-toupper(Methomlist$AGI_Number)
DDmetlist<-merge(DD,Methomlist,all.x=T,by.x="Gene",by.y="AGI_Number")

#################################################################
#add genes with high impact SNPs and allele frequency difference#
#################################################################
#SNPeffdfHIGH2<-merge(SNPeffdfHIGH,AFdata2)

HIGHlist<-SNPeffdfHIGH2[abs(SNPeffdfHIGH2$AF_Klet-SNPeffdfHIGH2$AF_Kowa)>0.4,]
Klet_UG<-read.table("Klet_UG.table",header=T,sep="\t")
Kowa_UG<-read.table("Kowa_UG.table",header=T,sep="\t")
KK_UG<-cbind(Klet_UG[,1:3],Kowa_UG[,3])
names(KK_UG)[4]<-"AC.1"

library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)

#load gene lists
Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)
Ahal_ortholist=read.delim("HalleriGeneOrthogroup.txt", header=TRUE)

# Load the file of Thaliana descriptions
Thal_description=read.delim("TAIR10_functional_descriptions.txt", header=TRUE)
colnames(Thal_description)[1]<-"Gene"
test1<-substr(Thal_description[,1],1,9)
Thal_description[,1]<-test1

## Load GFF files and add promotor region to genes

### Step 1: import genes from GFF files

lyr_txdm = makeTxDbFromGFF("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt", format = c("gff3"), organism = "Arabidopsis lyrata")
Ahal_txdm = makeTxDbFromGFF("annotation.gff", format = c("gff3"), organism = "Arabidopsis lyrata")

#View the information imported:
lyr_txdm
Ahal_txdm

#Define a function that allows the addition of promotor regions
extend <- function(x, upstream=0, downstream=0) #from https://support.bioconductor.org/p/78652/
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

# Use the function to extend promotor region 2kb
  # Upstream=2000 adds 2kb upstream from transcription start site
  # Downstream=0 adds 0bp downstream from the gene.
HalGenes_plusPromot = extend(genes(Ahal_txdm),upstream=2000, downstream=0)
LyrGenes_plusPromot = extend(genes(lyr_txdm),upstream=2000, downstream=0)

### Create a file with the list of Thaliana genes and descriptions

#lyrata genes to Orthogroup
HIGHlist_lyratagenes_df<-merge(HIGHlist,as.data.frame(LyrGenes_plusPromot),by.x="Gene",by.y="gene_id",all.x=T)
colnames(HIGHlist_lyratagenes_df)=c("Gene","CHROM", "POS", "Base","Characterization","Effect","AF_Klet","AF_Kowa","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand")
HIGHlist_lyratagenes_df$gene_start<-ifelse(HIGHlist_lyratagenes_df$gene_strand=="+",HIGHlist_lyratagenes_df$gene_start+2000,HIGHlist_lyratagenes_df$gene_start)
HIGHlist_lyratagenes_df$gene_end<-ifelse(HIGHlist_lyratagenes_df$gene_strand=="-",HIGHlist_lyratagenes_df$gene_end-2000,HIGHlist_lyratagenes_df$gene_end)
HIGHlist_lyratagenes_df$gene_size<-HIGHlist_lyratagenes_df$gene_size-2000


HIGHlist_lyratagenes_df_w_orthogroup= merge(HIGHlist_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)

#Orthogroup to Thaliana Genes
HIGHlist_lyratagenes_df_w_orthogroup_Thalgenes= merge(HIGHlist_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(HIGHlist_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(HIGHlist_lyratagenes_df_w_orthogroup_Thalgenes)[15]="Gene"

#Thaliana Genes to Gene Descriptions
Thalgenes_described_HIGHlist_lyratagenes_df_w_orthogroup= merge(Thal_description, HIGHlist_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)

Thalgenes_described_HIGHlist_lyratagenes_df_w_orthogroup$Gene<-toupper(Thalgenes_described_HIGHlist_lyratagenes_df_w_orthogroup$Gene)
Methomlist$AGI_Number<-toupper(Methomlist$AGI_Number)
Thalgenes_described_HIGHlist_lyratagenes_df_w_orthogroupmetlist<-merge(Thalgenes_described_HIGHlist_lyratagenes_df_w_orthogroup,Methomlist,all.x=T,by.x="Gene",by.y="AGI_Number")


#Generate file
write.table(Thalgenes_described_HIGHlist_lyratagenes_df_w_orthogroupmetlist,"GenesHIGHeffect_above40percentAFdiff_KK_arenosa.txt", sep="\t", row.names=F)

HighSNPeff<-read.table("GenesHIGHeffect_above40percentAFdiff_KK_arenosa.txt",sep="\t",header=T)
diff<-ifelse(paste(HighSNPeff$CHROM,HighSNPeff$POS,sep="_")%in%paste(KK_UG$CHROM,KK_UG$POS,sep="_"),"YES","NO")

summarydiff<-ifelse(diff=="YES",(KK_UG$AC[match(paste(HighSNPeff$CHROM,HighSNPeff$POS,sep="_"),paste(KK_UG$CHROM,KK_UG$POS,sep="_"))]/n1)-(KK_UG$AC.1[match(paste(HighSNPeff$CHROM,HighSNPeff$POS,sep="_"),paste(KK_UG$CHROM,KK_UG$POS,sep="_"))]/n2),"NO")
HighSNPeffdiff<-cbind(HighSNPeff,HighSNPeff$AF_Klet-HighSNPeff$AF_Kowa,summarydiff)
names(HighSNPeffdiff)[33:34]<-c("AF_Klet-AF_Kowa","UG_AF_diff")

TH_ortho_numbers<-table(Thal_ortholist$OrthoGroup)
TH_ortho_number<-data.frame(TH_ortho_numbers)
names(TH_ortho_number)<-c("OrthoGroup","Frequency")
Ah_ortho_numbers<-table(Ahal_ortholist$OrthoGroup)
Ah_ortho_number<-data.frame(Ah_ortho_numbers)
names(Ah_ortho_number)<-c("OrthoGroup","Frequency")
Ly_ortho_numbers<-table(Lyr_ortholist$OrthoGroup)
Ly_ortho_number<-data.frame(Ly_ortho_numbers)
names(Ly_ortho_number)<-c("OrthoGroup","Frequency")

Ah_paralogue_number<-Ah_ortho_number$Frequency[match(HighSNPeffdiff$OrthoGroup,Ah_ortho_number$OrthoGroup)]
Ath_paralogue_number<-TH_ortho_number$Frequency[match(HighSNPeffdiff$OrthoGroup,TH_ortho_number$OrthoGroup)]
Al_paralogue_number<-Ly_ortho_number$Frequency[match(HighSNPeffdiff$OrthoGroup,Ly_ortho_number$OrthoGroup)]

HighSNPeffdiff$Ah_paralogue_number<-Ah_paralogue_number
HighSNPeffdiff$Ath_paralogue_number<-Ath_paralogue_number
HighSNPeffdiff$Al_paralogue_number<-Al_paralogue_number


write.table(HighSNPeffdiff,"GenesHIGHeffect_above40percentAFdiff_withUG_KK_arenosa.txt", sep="\t", row.names=F)


options(java.parameters = "-Xmx24000m")

write.xlsx2(FSTmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="Fst",col.names=TRUE,row.names=FALSE)
write.xlsx2(DXYmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="Dxy",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(AFDabsmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="AFDabs",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(NIELSENmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="Nielsen",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(VARLDmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="VarLD",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(FLKmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="Flk",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(DDmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="DD",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(HighSNPeffdiff,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="Highimpact",col.names=TRUE,row.names=FALSE,append=T)

write.csv(FSTmetlist,"FST_KletKowa_01percent_arenosa.csv",row.names=F)
write.csv(DXYmetlist,"DXY_KletKowa_01percent_arenosa.csv",row.names=F)
write.csv(AFDabsmetlist,"AFDabs_KletKowa_01percent_arenosa.csv",row.names=F)
write.csv(NIELSENmetlist,"NIELSEN_KletKowa_01percent_arenosa.csv",row.names=F)
write.csv(VARLDmetlist,"VARLD_KletKowa_01percent_arenosa.csv",row.names=F)
write.csv(FLKmetlist,"FLK_KletKowa_01percent_arenosa.csv",row.names=F)
write.csv(DDmetlist,"DD_KletKowa_01percent_arenosa.csv",row.names=F)



FSTmetlist<-read.csv("FST_KletKowa_01percent_arenosa.csv",header=T)
DXYmetlist<-read.csv("DXY_KletKowa_01percent_arenosa.csv",header=T)
AFDabsmetlist<-read.csv("AFDabs_KletKowa_01percent_arenosa.csv",header=T)
NIELSENmetlist<-read.csv("NIELSEN_KletKowa_01percent_arenosa.csv",header=T)
VARLDmetlist<-read.csv("VARLD_KletKowa_01percent_arenosa.csv",header=T)
FLKmetlist<-read.csv("FLK_KletKowa_01percent_arenosa.csv",header=T)
DDmetlist<-read.csv("DD_KletKowa_01percent_arenosa.csv",header=T)
HighSNPeffdiff<-read.table("GenesHIGHeffect_above40percentAFdiff_withUG_KK_arenosa.txt", sep="\t", header=T)



