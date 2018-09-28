

###############
#Metrics plots#
###############

AFdata<-read.table("AllpostUGtests.csv",sep="\t",header=TRUE,encoding="UTF-8")
test2<-as.character(AFdata$Scaffold)
AFdata$Scaffold<-test2

require(xlsx)
Candidates2<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",1,header=TRUE)

Candidates3<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates<-Candidates3[!Candidates3$Gene=="",]

test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*1.5 #bp genemiddle +/- windowsize will be displayed


winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")
require(vegan)
TD_Klet$TajimasD_Klet<-decostand(TD_Klet$TajimasD_Klet,"standardize")
TD_Kowa$TajimasD_Kowa<-decostand(TD_Kowa$TajimasD_Kowa,"standardize")
TD_Klet$FayandWusH_Klet<-decostand(TD_Klet$FayandWusH_Klet,"standardize")
TD_Kowa$FayandWusH_Kowa<-decostand(TD_Kowa$FayandWusH_Kowa,"standardize")

test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

data1<-read.table("Klet_Sweed_scaffold1.table",sep="\t",header=T)
data1$Scaffold<-rep("Scaffold_1",nrow(data1))
data2<-read.table("Klet_Sweed_scaffold2.table",sep="\t",header=T)
data2$Scaffold<-rep("Scaffold_2",nrow(data2))
data3<-read.table("Klet_Sweed_scaffold3.table",sep="\t",header=T)
data3$Scaffold<-rep("Scaffold_3",nrow(data3))
data4<-read.table("Klet_Sweed_scaffold4.table",sep="\t",header=T)
data4$Scaffold<-rep("Scaffold_4",nrow(data4))
data5<-read.table("Klet_Sweed_scaffold5.table",sep="\t",header=T)
data5$Scaffold<-rep("Scaffold_5",nrow(data5))
data6<-read.table("Klet_Sweed_scaffold6.table",sep="\t",header=T)
data6$Scaffold<-rep("Scaffold_6",nrow(data6))
data7<-read.table("Klet_Sweed_scaffold7.table",sep="\t",header=T)
data7$Scaffold<-rep("Scaffold_7",nrow(data7))
data8<-read.table("Klet_Sweed_scaffold8.table",sep="\t",header=T)
data8$Scaffold<-rep("Scaffold_8",nrow(data8))

Sweed<-rbind(data1,data2,data3,data4,data5,data6,data7,data8)
Sweed<-Sweed[,c(4,1:3)]
Sweed[,1]<-as.factor(Sweed[,1])

test2<-as.character(tolower(Sweed$Scaffold))
Sweed$Scaffold<-test2
for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	interval<-which(AF1$Position>=(Candidates$gene_start[j]-windowsize[j])&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]))
	interval2<-c(min(interval)-1,interval,max(interval)+1)
	Sweedofinterest<-AF1[interval2,]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

dataS1<-read.table("Kowa_Sweed_scaffold1.table",sep="\t",header=T)
dataS1$Scaffold<-rep("Scaffold_1",nrow(dataS1))
dataS2<-read.table("Kowa_Sweed_scaffold2.table",sep="\t",header=T)
dataS2$Scaffold<-rep("Scaffold_2",nrow(dataS2))
dataS3<-read.table("Kowa_Sweed_scaffold3.table",sep="\t",header=T)
dataS3$Scaffold<-rep("Scaffold_3",nrow(dataS3))
dataS4<-read.table("Kowa_Sweed_scaffold4.table",sep="\t",header=T)
dataS4$Scaffold<-rep("Scaffold_4",nrow(dataS4))
dataS5<-read.table("Kowa_Sweed_scaffold5.table",sep="\t",header=T)
dataS5$Scaffold<-rep("Scaffold_5",nrow(dataS5))
dataS6<-read.table("Kowa_Sweed_scaffold6.table",sep="\t",header=T)
dataS6$Scaffold<-rep("Scaffold_6",nrow(dataS6))
dataS7<-read.table("Kowa_Sweed_scaffold7.table",sep="\t",header=T)
dataS7$Scaffold<-rep("Scaffold_7",nrow(dataS7))
dataS8<-read.table("Kowa_Sweed_scaffold8.table",sep="\t",header=T)
dataS8$Scaffold<-rep("Scaffold_8",nrow(dataS8))

Sweed_Kowa<-rbind(dataS1,dataS2,dataS3,dataS4,dataS5,dataS6,dataS7,dataS8)
Sweed_Kowa<-Sweed_Kowa[,c(4,1:3)]
Sweed_Kowa[,1]<-as.factor(Sweed_Kowa[,1])

test2<-as.character(tolower(Sweed_Kowa$Scaffold))
Sweed_Kowa$Scaffold<-test2
for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")


genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)


for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("FST_01percent_KletKowa_Metrics_arenosa_",j,"_",Candidates$Lyr_Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))

	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}



AFdata<-read.table("AllpostUGtests.csv",sep="\t",header=TRUE,encoding="UTF-8")
test2<-as.character(AFdata$Scaffold)
AFdata$Scaffold<-test2

Candidates2<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",4,header=TRUE)

Candidates3<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates<-Candidates3[!Candidates3$Gene=="",]
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)


test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*1.5 #bp genemiddle +/- windowsize will be displayed


winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")
require(vegan)
TD_Klet$TajimasD_Klet<-decostand(TD_Klet$TajimasD_Klet,"standardize")
TD_Kowa$TajimasD_Kowa<-decostand(TD_Kowa$TajimasD_Kowa,"standardize")
TD_Klet$FayandWusH_Klet<-decostand(TD_Klet$FayandWusH_Klet,"standardize")
TD_Kowa$FayandWusH_Kowa<-decostand(TD_Kowa$FayandWusH_Kowa,"standardize")


test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	Sweedofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")

genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)


for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("NIELSEN_01percent_KletKowa_Metrics_arenosa_",j,"_",Candidates$Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))

	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}


AFdata<-read.table("AllpostUGtests.csv",sep="\t",header=TRUE,encoding="UTF-8")
test2<-as.character(AFdata$Scaffold)
AFdata$Scaffold<-test2

Candidates2<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",2,header=TRUE)

Candidates3<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates<-Candidates3[!Candidates3$Gene=="",]

test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)


genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*1.5 #bp genemiddle +/- windowsize will be displayed

winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")
require(vegan)
TD_Klet$TajimasD_Klet<-decostand(TD_Klet$TajimasD_Klet,"standardize")
TD_Kowa$TajimasD_Kowa<-decostand(TD_Kowa$TajimasD_Kowa,"standardize")
TD_Klet$FayandWusH_Klet<-decostand(TD_Klet$FayandWusH_Klet,"standardize")
TD_Kowa$FayandWusH_Kowa<-decostand(TD_Kowa$FayandWusH_Kowa,"standardize")


test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	Sweedofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")

genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)


for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("DXY_01percent_KletKowa_Metrics_arenosa_",j,"_",Candidates$Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))


	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}



AFdata<-read.table("AllpostUGtests.csv",sep="\t",header=TRUE,encoding="UTF-8")
test2<-as.character(AFdata$Scaffold)
AFdata$Scaffold<-test2

Candidates2<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",6,header=TRUE)

Candidates3<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates<-Candidates3[!Candidates3$Gene=="",]

test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)


genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*1.5 #bp genemiddle +/- windowsize will be displayed

winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")
require(vegan)
TD_Klet$TajimasD_Klet<-decostand(TD_Klet$TajimasD_Klet,"standardize")
TD_Kowa$TajimasD_Kowa<-decostand(TD_Kowa$TajimasD_Kowa,"standardize")
TD_Klet$FayandWusH_Klet<-decostand(TD_Klet$FayandWusH_Klet,"standardize")
TD_Kowa$FayandWusH_Kowa<-decostand(TD_Kowa$FayandWusH_Kowa,"standardize")

test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	Sweedofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")
genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)



for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("FLK_01percent_KletKowa_Metrics_arenosa_",j,"_",Candidates$Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))

	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}



AFdata<-read.table("AllpostUGtests.csv",sep="\t",header=TRUE,encoding="UTF-8")
test2<-as.character(AFdata$Scaffold)
AFdata$Scaffold<-test2

Candidates2<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",5,header=TRUE)

Candidates3<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates<-Candidates3[!Candidates3$Gene=="",]

test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*1.5 #bp genemiddle +/- windowsize will be displayed

winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")
require(vegan)
TD_Klet$TajimasD_Klet<-decostand(TD_Klet$TajimasD_Klet,"standardize")
TD_Kowa$TajimasD_Kowa<-decostand(TD_Kowa$TajimasD_Kowa,"standardize")
TD_Klet$FayandWusH_Klet<-decostand(TD_Klet$FayandWusH_Klet,"standardize")
TD_Kowa$FayandWusH_Kowa<-decostand(TD_Kowa$FayandWusH_Kowa,"standardize")


test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	Sweedofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")
genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)



for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("VARLD_01percent_KletKowa_Metrics_arenosa_",j,"_",Candidates$Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))


	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}


Candidates2<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",3,header=TRUE)

Candidates3<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates<-Candidates3[!Candidates3$Gene=="",]

test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)


genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*1.5 #bp genemiddle +/- windowsize will be displayed

winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")
require(vegan)
TD_Klet$TajimasD_Klet<-decostand(TD_Klet$TajimasD_Klet,"standardize")
TD_Kowa$TajimasD_Kowa<-decostand(TD_Kowa$TajimasD_Kowa,"standardize")
TD_Klet$FayandWusH_Klet<-decostand(TD_Klet$FayandWusH_Klet,"standardize")
TD_Kowa$FayandWusH_Kowa<-decostand(TD_Kowa$FayandWusH_Kowa,"standardize")


test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	Sweedofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")

genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)


for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("AFDabs_01percent_KletKowa_Metrics_arenosa_",j,"_",Candidates$Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))


	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}


Candidatesall<-read.csv("Lyratagenes_01percent_KletKowa_reflyrata_arenosa.csv",header=TRUE,sep=";")

Candidates<-Candidatesall

Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)


test2<-as.character(Candidates$Lyr_Gene)
Candidates$Lyr_Gene<-test2
genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*1.5 #bp genemiddle +/- windowsize will be displayed



winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")
require(vegan)
TD_Klet$TajimasD_Klet<-decostand(TD_Klet$TajimasD_Klet,"standardize")
TD_Kowa$TajimasD_Kowa<-decostand(TD_Kowa$TajimasD_Kowa,"standardize")
TD_Klet$FayandWusH_Klet<-decostand(TD_Klet$FayandWusH_Klet,"standardize")
TD_Kowa$FayandWusH_Kowa<-decostand(TD_Kowa$FayandWusH_Kowa,"standardize")


test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	Sweedofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")

genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)


for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("Lyrataonly_01percent_KletKowa_Metrics_arenosa_",j,"_",Candidates$Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))


	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}


Candidates2<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",7,header=TRUE)

Candidates3<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates<-Candidates3[!Candidates3$Gene=="",]

test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*1.5 #bp genemiddle +/- windowsize will be displayed

winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")

test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	Sweedofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")
genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)



for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("DD_01percent_KletKowa_Metrics_",j,"_",Candidates$Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))

	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}





AFdata<-read.table("AllpostUGtests.csv",sep="\t",header=TRUE,encoding="UTF-8")
names(AFdata)<-c("Scaffold","Start_pos","End_pos","Fst","DD","Nielsen","Dxy","Flk","Flk_p_value","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AFKrom","AFKosi","TajimasD_Krom","TajimasD_Kosi","FayandWusH_Krom","FayandWusH_Kosi")
test2<-as.character(AFdata$Scaffold)
AFdata$Scaffold<-test2

Candidates2<-read.xlsx2("Genes_01percent_KletKowa_arenosa.xlsx",3,header=TRUE)

Candidates3<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates<-Candidates3[!Candidates3$Gene=="",]

test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*4 #bp genemiddle +/- windowsize will be displayed
test3<-as.character(Candidates$Chr)
Candidates$Chr<-test3

winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")

test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	Sweedofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")

genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)


for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("AFDabs_01percent_KletKowa_Metrics_",j,"_",Candidates$Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))

	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}







AFdata<-read.table("AllpostUGtests.csv",sep="\t",header=TRUE,encoding="UTF-8")
names(AFdata)<-c("Scaffold","Start_pos","End_pos","Fst","DD","Nielsen","Dxy","Flk","Flk_p_value","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AFKrom","AFKosi","TajimasD_Krom","TajimasD_Kosi","FayandWusH_Krom","FayandWusH_Kosi")
test2<-as.character(AFdata$Scaffold)
AFdata$Scaffold<-test2

Candidatesall<-read.csv("Summary_genes_1percent_KletKowa_TD.csv",header=TRUE,sep=";")
Candidates2<-Candidatesall[Candidatesall$TajimasD_Klet<=-3.24163863346625,]
Candidates3<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates<-Candidates3[!Candidates3$Gene=="",]

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*4 #bp genemiddle +/- windowsize will be displayed
test3<-as.character(Candidates$Chr)
Candidates$Chr<-test3

winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")

test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	Sweedofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")

genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)


for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("TD_01percent_KletKowa_Metrics_",j,"_",Candidates$Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))


	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}





##########################################################################
#Candidate genes all metrics 0.1% lyrata genes without thaliana orthologues#
##########################################################################
AFdata<-read.table("AllpostUGtests.csv",sep="\t",header=TRUE,encoding="UTF-8")
names(AFdata)<-c("Scaffold","Start_pos","End_pos","Fst","DD","Nielsen","Dxy","Flk","Flk_p_value","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AFKrom","AFKosi","TajimasD_Krom","TajimasD_Kosi","FayandWusH_Krom","FayandWusH_Kosi")
test2<-as.character(AFdata$Scaffold)
AFdata$Scaffold<-test2

Candidates<-read.csv("Lyratagenes_01percent_KletKowa_reflyrata.csv",header=TRUE,sep=";")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*4 #bp genemiddle +/- windowsize will be displayed
test3<-as.character(Candidates$Chr)
Candidates$Chr<-test3

winmiddle<-c((AFdata$Start_pos+AFdata$End_pos)/2)
AFdata<-cbind(AFdata,winmiddle)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Scaffold==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")

TD_Klet<-read.csv("TDFWH_Klet.csv",header=TRUE,sep="\t")
TD_Kowa<-read.csv("TDFWH_Kowa.csv",header=TRUE,sep="\t")

test2<-as.character(TD_Klet$Scaffold)
TD_Klet$Scaffold<-test2
winmiddleTDKlet<-c((TD_Klet$Start_pos+TD_Klet$End_pos)/2)
TD_Klet<-cbind(TD_Klet,winmiddleTDKlet)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKletofinterest", j, sep ="")
 	AF1<-TD_Klet[TD_Klet$Scaffold==Candidates$Chr[j],]
	TDKletofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKletofinterest)
	}
rm(TDKletofinterest)
list<-ls(pattern="^TDKletofinterest")

test2<-as.character(TD_Kowa$Scaffold)
TD_Kowa$Scaffold<-test2
winmiddleTDKowa<-c((TD_Kowa$Start_pos+TD_Kowa$End_pos)/2)
TD_Kowa<-cbind(TD_Kowa,winmiddleTDKowa)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKowaofinterest", j, sep ="")
 	AF1<-TD_Kowa[TD_Kowa$Scaffold==Candidates$Chr[j],]
	TDKowaofinterest<-AF1[AF1$winmiddle>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$winmiddle<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKowaofinterest)
	}
rm(TDKowaofinterest)
list<-ls(pattern="^TDKowaofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_perc<-quantile(AFdata$DD,0.001)
Fst_perc<-quantile(AFdata$Fst,0.999)
Nielsen_perc<-quantile(AFdata$Nielsen,0.999)
Dxy_perc<-quantile(AFdata$Dxy,0.999)
Flk_perc<-quantile(AFdata$Flk,0.999)
VarLD_perc<-quantile(AFdata$VarLD,0.999)
TDKlet_perc<-quantile(TD_Klet$TajimasD_Klet,0.001)
TDKowa_perc<-quantile(TD_Kowa$TajimasD_Kowa,0.001)

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

data1<-read.table("Klet_Sweed_scaffold1.table",sep="\t",header=T)
data1$Scaffold<-rep("Scaffold_1",nrow(data1))
data2<-read.table("Klet_Sweed_scaffold2.table",sep="\t",header=T)
data2$Scaffold<-rep("Scaffold_2",nrow(data2))
data3<-read.table("Klet_Sweed_scaffold3.table",sep="\t",header=T)
data3$Scaffold<-rep("Scaffold_3",nrow(data3))
data4<-read.table("Klet_Sweed_scaffold4.table",sep="\t",header=T)
data4$Scaffold<-rep("Scaffold_4",nrow(data4))
data5<-read.table("Klet_Sweed_scaffold5.table",sep="\t",header=T)
data5$Scaffold<-rep("Scaffold_5",nrow(data5))
data6<-read.table("Klet_Sweed_scaffold6.table",sep="\t",header=T)
data6$Scaffold<-rep("Scaffold_6",nrow(data6))
data7<-read.table("Klet_Sweed_scaffold7.table",sep="\t",header=T)
data7$Scaffold<-rep("Scaffold_7",nrow(data7))
data8<-read.table("Klet_Sweed_scaffold8.table",sep="\t",header=T)
data8$Scaffold<-rep("Scaffold_8",nrow(data8))

Sweed<-rbind(data1,data2,data3,data4,data5,data6,data7,data8)
Sweed<-Sweed[,c(4,1:3)]
Sweed[,1]<-as.factor(Sweed[,1])

test2<-as.character(tolower(Sweed$Scaffold))
Sweed$Scaffold<-test2
for (j in 1:nrow(Candidates))
	{nam <- paste("Sweedofinterest", j, sep ="")
 	AF1<-Sweed[Sweed$Scaffold==Candidates$Chr[j],]
	Sweedofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweedofinterest)
	}
rm(Sweedofinterest)
list<-ls(pattern="^Sweedofinterest")

dataS1<-read.table("Kowa_Sweed_scaffold1.table",sep="\t",header=T)
dataS1$Scaffold<-rep("Scaffold_1",nrow(dataS1))
dataS2<-read.table("Kowa_Sweed_scaffold2.table",sep="\t",header=T)
dataS2$Scaffold<-rep("Scaffold_2",nrow(dataS2))
dataS3<-read.table("Kowa_Sweed_scaffold3.table",sep="\t",header=T)
dataS3$Scaffold<-rep("Scaffold_3",nrow(dataS3))
dataS4<-read.table("Kowa_Sweed_scaffold4.table",sep="\t",header=T)
dataS4$Scaffold<-rep("Scaffold_4",nrow(dataS4))
dataS5<-read.table("Kowa_Sweed_scaffold5.table",sep="\t",header=T)
dataS5$Scaffold<-rep("Scaffold_5",nrow(dataS5))
dataS6<-read.table("Kowa_Sweed_scaffold6.table",sep="\t",header=T)
dataS6$Scaffold<-rep("Scaffold_6",nrow(dataS6))
dataS7<-read.table("Kowa_Sweed_scaffold7.table",sep="\t",header=T)
dataS7$Scaffold<-rep("Scaffold_7",nrow(dataS7))
dataS8<-read.table("Kowa_Sweed_scaffold8.table",sep="\t",header=T)
dataS8$Scaffold<-rep("Scaffold_8",nrow(dataS8))

Sweed_Kowa<-rbind(dataS1,dataS2,dataS3,dataS4,dataS5,dataS6,dataS7,dataS8)
Sweed_Kowa<-Sweed_Kowa[,c(4,1:3)]
Sweed_Kowa[,1]<-as.factor(Sweed_Kowa[,1])

test2<-as.character(tolower(Sweed_Kowa$Scaffold))
Sweed_Kowa$Scaffold<-test2
for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kowaofinterest", j, sep ="")
 	AF1<-Sweed_Kowa[Sweed_Kowa$Scaffold==Candidates$Chr[j],]
	Sweed_Kowaofinterest<-AF1[AF1$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kowaofinterest)
	}
rm(Sweed_Kowaofinterest)
list<-ls(pattern="^Sweed_Kowaofinterest")

genelist2<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist2)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
exonlist<-genelist2[genelist2$Type=="exon",]
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
IDlist<-strsplit(as.character(exonlist[,9]),"=")
IDs<-matrix(unlist(IDlist),ncol=3,byrow=TRUE)
exonlist$ID<-IDs[,3]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Lyr_Gene<-as.character(Candidates$Lyr_Gene)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Lyr_Gene[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)



for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("Lyrataonly01percent_KletKowa_Metrics_",j,"_",Candidates$Lyr_Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=20, height=26, units="cm", res=1000)
	par(mfrow=c(9,1))
	par(mar=c(1,5,1,1)+0.1)
	par(oma=c(4,2,0,0))

	plot(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(DD~residuals)),ylim=c(min(AFdata$DD),max(AFdata$DD)+0.15),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple")
	rect(Candidates$gene_start[j],min(AFdata$DD-0.3),Candidates$gene_end[j],max(AFdata$DD)+0.4,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$DD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="purple")
	abline(h=DD_abs_perc,lty=2,col="purple")
	text(x=genemiddle[[1]][j],y=max(AFdata$DD)+0.12,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=3)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdata$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=3)

	plot(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Fst)),ylim=c(min(AFdata$Fst),max(AFdata$Fst)+0.6),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink")
	rect(Candidates$gene_start[j],min(AFdata$Fst)-0.7,Candidates$gene_end[j],max(AFdata$Fst)+0.7,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Fst~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="hotpink")
	abline(h=Fst_perc,lty=2,col="hotpink")
	text(x=genemiddle[[1]][j],y=max(AFdata$Fst)+0.55,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=3)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdata$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=3)

	plot(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Nielsen)),ylim=c(min(AFdata$Nielsen),max(AFdata$Nielsen)+150),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red")
	rect(Candidates$gene_start[j],min(AFdata$Nielsen)-200,Candidates$gene_end[j],max(AFdata$Nielsen)+200,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Nielsen~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="red")
	abline(h=Nielsen_perc,lty=2,col="red")
	text(x=genemiddle[[1]][j],y=max(AFdata$Nielsen)+125,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=3)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdata$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=3)

	plot(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Dxy)),ylim=c(min(AFdata$Dxy),max(AFdata$Dxy)+0.8),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="orange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="orange")
	rect(Candidates$gene_start[j],min(AFdata$Dxy-0.9),Candidates$gene_end[j],max(AFdata$Dxy)+0.9,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Dxy~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="orange")
	abline(h=Dxy_perc,lty=2,col="orange")
	text(x=genemiddle[[1]][j],y=max(AFdata$Dxy)+0.7,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=3)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdata$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=3)

	plot(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(Flk)),ylim=c(min(AFdata$Flk),max(AFdata$Flk)+5),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="yellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow")
	rect(Candidates$gene_start[j],min(AFdata$Flk)-5,Candidates$gene_end[j],max(AFdata$Flk)+8,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$Flk~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="yellow")
	abline(h=Flk_perc,lty=2,col="yellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$Flk)+4,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=3)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdata$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=3)

	plot(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,ylab=expression(bold(VarLD)),ylim=c(min(AFdata$VarLD),max(AFdata$VarLD)+50),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="greenyellow",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="greenyellow")
	rect(Candidates$gene_start[j],min(AFdata$VarLD)-60,Candidates$gene_end[j],max(AFdata$VarLD)+60,col="grey",border = NA)	
	lines(get(paste("AFofinterest",j,sep=""))$VarLD~get(paste("AFofinterest",j,sep=""))$winmiddle,lwd=2,col="greenyellow")
	abline(h=VarLD_perc,lty=2,col="greenyellow")
	text(x=genemiddle[[1]][j],y=max(AFdata$VarLD)+45,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdata$VarLD)+30,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdata$VarLD)+30,code=arrowdir[i],length=0.1,col="grey",lwd=3)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdata$VarLD)+30,x1=Candidates$gene_end[j],y1=max(AFdata$VarLD)+30,code=arrowdir2[j],length=0.1,col="red",lwd=3)

	plot(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(TD~Klet)),ylim=c(min(TD_Klet$TajimasD_Klet),max(TD_Klet$TajimasD_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen")
	rect(Candidates$gene_start[j],min(TD_Klet$TajimasD_Klet)-60,Candidates$gene_end[j],max(TD_Klet$TajimasD_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$TajimasD_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkgreen")
	lines(get(paste("TDKowaofinterest",j,sep=""))$TajimasD_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightgreen")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$TajimasD_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$TajimasD_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=3)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$TajimasD_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$TajimasD_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=3)
	mtext(text=expression(bold(TD~Kowa)),line=5,side=2,col="lightgreen",cex=1)

	plot(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,ylab=expression(bold(FWH~Klet)),ylim=c(min(TD_Klet$FayandWusH_Klet),max(TD_Klet$FayandWusH_Klet)+3),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="darkblue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkblue")
	rect(Candidates$gene_start[j],min(TD_Klet$FayandWusH_Klet)-60,Candidates$gene_end[j],max(TD_Klet$FayandWusH_Klet)+60,col="grey",border = NA)	
	lines(get(paste("TDKletofinterest",j,sep=""))$FayandWusH_Klet~get(paste("TDKletofinterest",j,sep=""))$winmiddle,lwd=2,col="darkblue")
	lines(get(paste("TDKowaofinterest",j,sep=""))$FayandWusH_Kowa~get(paste("TDKowaofinterest",j,sep=""))$winmiddle,lwd=2,col="lightblue")
	text(x=genemiddle[[1]][j],y=max(TD_Klet$FayandWusH_Klet)+2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir[i],length=0.1,col="grey",lwd=3)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Klet$FayandWusH_Klet)+1,x1=Candidates$gene_end[j],y1=max(TD_Klet$FayandWusH_Klet)+1,code=arrowdir2[j],length=0.1,col="red",lwd=3)
	mtext(text=expression(bold(FWH~Kowa)),line=5,side=2,col="lightblue",cex=0.99)

	plot(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,ylab=expression(bold(Sweed~Klet)),ylim=c(min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4")
	rect(Candidates$gene_start[j],min(Sweed$Likelihood)-60,Candidates$gene_end[j],max(Sweed$Likelihood)+60,col="grey",border = NA)	
	lines(get(paste("Sweedofinterest",j,sep=""))$Likelihood~get(paste("Sweedofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kowaofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
	text(x=genemiddle[[1]][j],y=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=3)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("Sweedofinterest",j,sep=""))$Likelihood),max(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweedofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Kowaofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=3)
	mtext(text=expression(bold(Sweed~Kowa)),line=5,side=2,col="lightsalmon",cex=0.99)
	mtext(text=expression(bold(Scaffold~position~(bp))),side=1,line=3,outer=TRUE,cex=1.5)


	dev.off()
	}


