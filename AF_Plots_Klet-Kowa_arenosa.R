
####################################
#Folded AF plots with Klet and Kowa#
####################################
n1<-36 #number of Klet individuals * ploidy
n2<-32 #number of Kowa individuals * ploidy


AFdata<-read.table("LyrataKletKowaGS_woGQfil_arenosa.csv",header=TRUE,sep="\t")

test<-as.character(AFdata$Chrom)
AFdata[,1]<-test

AFdata[,3]<-AFdata[,3]/n1
AFdata[,6]<-AFdata[,6]/n2

##########################
#Candidate genes Fst 0.1%#
##########################
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

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Chrom==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$POS>=(genemiddle[[1]][j]-windowsize[j])&AF1$POS<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")


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
rm(Genesofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

test<-as.factor(AFdata$Chrom)
AFdata[,1]<-test

winsize=25
AFmatrix=data.frame()
for (i in levels(AFdata[,1]))
        {datanow<-AFdata[AFdata[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,6])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix=rbind(AFmatrix,wmatrix)
        }
names(AFmatrix)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Klet","AF_Kowa")


test<-as.character(AFdata$Chrom)
AFdata[,1]<-test
test<-as.character(AFmatrix$Chrom)
AFmatrix[,1]<-test

test2<-AFmatrix
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix[,i]))
	}
AFmatrix<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest", j, sep ="")
 	AF1<-AFmatrix[AFmatrix$Chrom==Candidates$Chr[j],]
	AFfreqofinterest<-AF1[AF1$Window_start>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_start<=(genemiddle[[1]][j]+windowsize[j]+100000)|AF1$Window_end>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_end<=(genemiddle[[1]][j]+windowsize[j]+100000),]
	assign(nam, AFfreqofinterest)
	}
rm(AFfreqofinterest)

SNPeffFst<-read.table("SNPeff_1perc_KK_arenosa.table",header=T,sep="\t")

AFdatasnp<-cbind(AFdata[,1:3],AFdata[,6])
names(AFdatasnp)<-c("CHROM","POS","AC","AC.1")
SNPeffdfAF<-merge(SNPeffFst,AFdatasnp)
SNPeffdfAFo<-SNPeffdfAF[SNPeffdfAF$Effect=="MODERATE",]
SNPeffdfAFp<-SNPeffdfAF[SNPeffdfAF$Effect=="MODIFIER",]
SNPeffdfAFr<-SNPeffdfAF[SNPeffdfAF$Effect=="HIGH",]

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
	nam <- paste("Fst01percent_AFKletKowa_absoluteAFplots_candidate_genes_arenosa_",j,"_",Candidates$Lyr_Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=18, height=18, units="cm", res=1000)
	plot(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS,ylab="AF difference Klet-Kowa",xlab="Scaffold position (bp)",ylim=c(0,1.4),cex.lab=1.5,xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]))
	rect(Candidates$gene_start[j],-2, Candidates$gene_end[j], 2,col="grey",border = NA)	
	points(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS)
	text(x=genemiddle[[1]][j],y=1.2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=1.1,x1=Candidates$gene_end[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Klet[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="darkblue")
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Kowa[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="lightblue")
	legend("top",legend=c("25SNP-window AF Klet","25SNP-window AF Kowa"),pch=15,col=c("darkblue","lightblue"),text.font=2,bg="white")
	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])],col="purple")
	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])],col="orange")
	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])],col="red")
	legend("topleft",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2)
	dev.off()
	}


#########################
#Candidate genes DD 0.1%#
#########################
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

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Chrom==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$POS>=(genemiddle[[1]][j]-windowsize[j])&AF1$POS<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")


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
rm(Genesofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

test<-as.factor(AFdata$Chrom)
AFdata[,1]<-test

winsize=25
AFmatrix=data.frame()
for (i in levels(AFdata[,1]))
        {datanow<-AFdata[AFdata[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,6])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix=rbind(AFmatrix,wmatrix)
        }
names(AFmatrix)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Klet","AF_Kowa")


test<-as.character(AFdata$Chrom)
AFdata[,1]<-test
test<-as.character(AFmatrix$Chrom)
AFmatrix[,1]<-test

test2<-AFmatrix
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix[,i]))
	}
AFmatrix<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest", j, sep ="")
 	AF1<-AFmatrix[AFmatrix$Chrom==Candidates$Chr[j],]
	AFfreqofinterest<-AF1[AF1$Window_start>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_start<=(genemiddle[[1]][j]+windowsize[j]+100000)|AF1$Window_end>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_end<=(genemiddle[[1]][j]+windowsize[j]+100000),]
	assign(nam, AFfreqofinterest)
	}
rm(AFfreqofinterest)

SNPeffDD<-SNPeffFst

AFdatasnp<-cbind(AFdata[,1:3],AFdata[,6])
names(AFdatasnp)<-c("CHROM","POS","AC","AC.1")
SNPeffdfAF<-merge(SNPeffDD,AFdatasnp)
SNPeffdfAFo<-SNPeffdfAF[SNPeffdfAF$Effect=="MODERATE",]
SNPeffdfAFp<-SNPeffdfAF[SNPeffdfAF$Effect=="MODIFIER",]
SNPeffdfAFr<-SNPeffdfAF[SNPeffdfAF$Effect=="HIGH",]
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
	nam <- paste("DD01percent_AFKletKowa_absoluteAFplots_candidate_genes_arenosa_",j,"_",Candidates$Lyr_Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=18, height=18, units="cm", res=1000)
	plot(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS,ylab="AF difference Klet-Kowa",xlab="Scaffold position (bp)",ylim=c(0,1.4),cex.lab=1.5,xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]))
	rect(Candidates$gene_start[j],-2, Candidates$gene_end[j], 2,col="grey",border = NA)	
	points(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS)
	text(x=genemiddle[[1]][j],y=1.2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=1.1,x1=Candidates$gene_end[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Klet[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="darkblue")
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Kowa[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="lightblue")
	legend("top",legend=c("25SNP-window AF Klet","25SNP-window AF Kowa"),pch=15,col=c("darkblue","lightblue"),text.font=2,bg="white")
	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])],col="purple")
	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])],col="orange")
	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])],col="red")
	legend("topleft",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2)
	dev.off()
	}

##############################
#Candidate genes Nielsen 0.1%#
##############################

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

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Chrom==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$POS>=(genemiddle[[1]][j]-windowsize[j])&AF1$POS<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")


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
rm(Genesofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

test<-as.factor(AFdata$Chrom)
AFdata[,1]<-test

winsize=25
AFmatrix=data.frame()
for (i in levels(AFdata[,1]))
        {datanow<-AFdata[AFdata[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,6])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix=rbind(AFmatrix,wmatrix)
        }
names(AFmatrix)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Klet","AF_Kowa")


test<-as.character(AFdata$Chrom)
AFdata[,1]<-test
test<-as.character(AFmatrix$Chrom)
AFmatrix[,1]<-test

test2<-AFmatrix
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix[,i]))
	}
AFmatrix<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest", j, sep ="")
 	AF1<-AFmatrix[AFmatrix$Chrom==Candidates$Chr[j],]
	AFfreqofinterest<-AF1[AF1$Window_start>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_start<=(genemiddle[[1]][j]+windowsize[j]+100000)|AF1$Window_end>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_end<=(genemiddle[[1]][j]+windowsize[j]+100000),]
	assign(nam, AFfreqofinterest)
	}
rm(AFfreqofinterest)

SNPeffNielsen<-SNPeffFst

AFdatasnp<-cbind(AFdata[,1:3],AFdata[,6])
names(AFdatasnp)<-c("CHROM","POS","AC","AC.1")
SNPeffdfAF<-merge(SNPeffNielsen,AFdatasnp)
SNPeffdfAFo<-SNPeffdfAF[SNPeffdfAF$Effect=="MODERATE",]
SNPeffdfAFp<-SNPeffdfAF[SNPeffdfAF$Effect=="MODIFIER",]
SNPeffdfAFr<-SNPeffdfAF[SNPeffdfAF$Effect=="HIGH",]
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
	nam <- paste("Nielsen01percent_AFKletKowa_absoluteAFplots_candidate_genes_arenosa_",j,"_",Candidates$Lyr_Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=18, height=18, units="cm", res=1000)
	plot(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS,ylab="AF difference Klet-Kowa",xlab="Scaffold position (bp)",ylim=c(0,1.4),cex.lab=1.5,xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]))
	rect(Candidates$gene_start[j],-2, Candidates$gene_end[j], 2,col="grey",border = NA)	
	points(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS)
	text(x=genemiddle[[1]][j],y=1.2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=1.1,x1=Candidates$gene_end[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Klet[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="darkblue")
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Kowa[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="lightblue")
	legend("top",legend=c("25SNP-window AF Klet","25SNP-window AF Kowa"),pch=15,col=c("darkblue","lightblue"),text.font=2,bg="white")
	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])],col="purple")
	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])],col="orange")
	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])],col="red")
	legend("topleft",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2)
	dev.off()
	}
##########################
#Candidate genes Dxy 0.1%#
##########################
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

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Chrom==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$POS>=(genemiddle[[1]][j]-windowsize[j])&AF1$POS<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")


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
rm(Genesofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

test<-as.factor(AFdata$Chrom)
AFdata[,1]<-test

winsize=25
AFmatrix=data.frame()
for (i in levels(AFdata[,1]))
        {datanow<-AFdata[AFdata[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,6])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix=rbind(AFmatrix,wmatrix)
        }
names(AFmatrix)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Klet","AF_Kowa")


test<-as.character(AFdata$Chrom)
AFdata[,1]<-test
test<-as.character(AFmatrix$Chrom)
AFmatrix[,1]<-test

test2<-AFmatrix
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix[,i]))
	}
AFmatrix<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest", j, sep ="")
 	AF1<-AFmatrix[AFmatrix$Chrom==Candidates$Chr[j],]
	AFfreqofinterest<-AF1[AF1$Window_start>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_start<=(genemiddle[[1]][j]+windowsize[j]+100000)|AF1$Window_end>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_end<=(genemiddle[[1]][j]+windowsize[j]+100000),]
	assign(nam, AFfreqofinterest)
	}
rm(AFfreqofinterest)

SNPeffDxy<-SNPeffFst

AFdatasnp<-cbind(AFdata[,1:3],AFdata[,6])
names(AFdatasnp)<-c("CHROM","POS","AC","AC.1")
SNPeffdfAF<-merge(SNPeffDxy,AFdatasnp)
SNPeffdfAFo<-SNPeffdfAF[SNPeffdfAF$Effect=="MODERATE",]
SNPeffdfAFp<-SNPeffdfAF[SNPeffdfAF$Effect=="MODIFIER",]
SNPeffdfAFr<-SNPeffdfAF[SNPeffdfAF$Effect=="HIGH",]
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
	nam <- paste("Dxy01percent_AFKletKowa_absoluteAFplots_candidate_genes_arenosa_",j,"_",Candidates$Lyr_Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=18, height=18, units="cm", res=1000)
	plot(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS,ylab="AF difference Klet-Kowa",xlab="Scaffold position (bp)",ylim=c(0,1.4),cex.lab=1.5,xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]))
	rect(Candidates$gene_start[j],-2, Candidates$gene_end[j], 2,col="grey",border = NA)	
	points(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS)
	text(x=genemiddle[[1]][j],y=1.2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=1.1,x1=Candidates$gene_end[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Klet[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="darkblue")
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Kowa[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="lightblue")
	legend("top",legend=c("25SNP-window AF Klet","25SNP-window AF Kowa"),pch=15,col=c("darkblue","lightblue"),text.font=2,bg="white")
	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])],col="purple")
	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])],col="orange")
	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])],col="red")
	legend("topleft",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2)
	dev.off()
	}

###########################
#Candidate genes Flk <0.05#
###########################
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

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Chrom==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$POS>=(genemiddle[[1]][j]-windowsize[j])&AF1$POS<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")


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
rm(Genesofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

test<-as.factor(AFdata$Chrom)
AFdata[,1]<-test

winsize=25
AFmatrix=data.frame()
for (i in levels(AFdata[,1]))
        {datanow<-AFdata[AFdata[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,6])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix=rbind(AFmatrix,wmatrix)
        }
names(AFmatrix)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Klet","AF_Kowa")


test<-as.character(AFdata$Chrom)
AFdata[,1]<-test
test<-as.character(AFmatrix$Chrom)
AFmatrix[,1]<-test

test2<-AFmatrix
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix[,i]))
	}
AFmatrix<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest", j, sep ="")
 	AF1<-AFmatrix[AFmatrix$Chrom==Candidates$Chr[j],]
	AFfreqofinterest<-AF1[AF1$Window_start>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_start<=(genemiddle[[1]][j]+windowsize[j]+100000)|AF1$Window_end>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_end<=(genemiddle[[1]][j]+windowsize[j]+100000),]
	assign(nam, AFfreqofinterest)
	}
rm(AFfreqofinterest)

SNPeffFlk<-SNPeffFst

AFdatasnp<-cbind(AFdata[,1:3],AFdata[,6])
names(AFdatasnp)<-c("CHROM","POS","AC","AC.1")
SNPeffdfAF<-merge(SNPeffFlk,AFdatasnp)
SNPeffdfAFo<-SNPeffdfAF[SNPeffdfAF$Effect=="MODERATE",]
SNPeffdfAFp<-SNPeffdfAF[SNPeffdfAF$Effect=="MODIFIER",]
SNPeffdfAFr<-SNPeffdfAF[SNPeffdfAF$Effect=="HIGH",]
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
	nam <- paste("Flk01percent_AFKletKowa_absoluteAFplots_candidate_genes_arenosa_",j,"_",Candidates$Lyr_Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=18, height=18, units="cm", res=1000)
	plot(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS,ylab="AF difference Klet-Kowa",xlab="Scaffold position (bp)",ylim=c(0,1.4),cex.lab=1.5,xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]))
	rect(Candidates$gene_start[j],-2, Candidates$gene_end[j], 2,col="grey",border = NA)	
	points(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS)
	text(x=genemiddle[[1]][j],y=1.2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=1.1,x1=Candidates$gene_end[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Klet[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="darkblue")
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Kowa[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="lightblue")
	legend("top",legend=c("25SNP-window AF Klet","25SNP-window AF Kowa"),pch=15,col=c("darkblue","lightblue"),text.font=2,bg="white")
	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])],col="purple")
	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])],col="orange")
	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])],col="red")
	legend("topleft",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2)
	dev.off()
	}






############################
#Candidate genes VarLD 0.1%#
############################
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

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Chrom==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$POS>=(genemiddle[[1]][j]-windowsize[j])&AF1$POS<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")


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
rm(Genesofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

test<-as.factor(AFdata$Chrom)
AFdata[,1]<-test

winsize=25
AFmatrix=data.frame()
for (i in levels(AFdata[,1]))
        {datanow<-AFdata[AFdata[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,6])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix=rbind(AFmatrix,wmatrix)
        }
names(AFmatrix)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Klet","AF_Kowa")


test<-as.character(AFdata$Chrom)
AFdata[,1]<-test
test<-as.character(AFmatrix$Chrom)
AFmatrix[,1]<-test

test2<-AFmatrix
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix[,i]))
	}
AFmatrix<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest", j, sep ="")
 	AF1<-AFmatrix[AFmatrix$Chrom==Candidates$Chr[j],]
	AFfreqofinterest<-AF1[AF1$Window_start>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_start<=(genemiddle[[1]][j]+windowsize[j]+100000)|AF1$Window_end>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_end<=(genemiddle[[1]][j]+windowsize[j]+100000),]
	assign(nam, AFfreqofinterest)
	}
rm(AFfreqofinterest)

SNPeffVarLD<-SNPeffFst

AFdatasnp<-cbind(AFdata[,1:3],AFdata[,6])
names(AFdatasnp)<-c("CHROM","POS","AC","AC.1")
SNPeffdfAF<-merge(SNPeffVarLD,AFdatasnp)
SNPeffdfAFo<-SNPeffdfAF[SNPeffdfAF$Effect=="MODERATE",]
SNPeffdfAFp<-SNPeffdfAF[SNPeffdfAF$Effect=="MODIFIER",]
SNPeffdfAFr<-SNPeffdfAF[SNPeffdfAF$Effect=="HIGH",]
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
	nam <- paste("VarLD01percent_AFKletKowa_absoluteAFplots_candidate_genes_arenosa_",j,"_",Candidates$Lyr_Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=18, height=18, units="cm", res=1000)
	plot(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS,ylab="AF difference Klet-Kowa",xlab="Scaffold position (bp)",ylim=c(0,1.4),cex.lab=1.5,xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]))
	rect(Candidates$gene_start[j],-2, Candidates$gene_end[j], 2,col="grey",border = NA)	
	points(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS)
	text(x=genemiddle[[1]][j],y=1.2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=1.1,x1=Candidates$gene_end[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Klet[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="darkblue")
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Kowa[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="lightblue")
	legend("top",legend=c("25SNP-window AF Klet","25SNP-window AF Kowa"),pch=15,col=c("darkblue","lightblue"),text.font=2,bg="white")
	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])],col="purple")
	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])],col="orange")
	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])],col="red")
	legend("topleft",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2)
	dev.off()
	}
#############################
#Candidate genes AFDabs 0.1%#
#############################
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

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Chrom==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$POS>=(genemiddle[[1]][j]-windowsize[j])&AF1$POS<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")


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
rm(Genesofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

test<-as.factor(AFdata$Chrom)
AFdata[,1]<-test

winsize=25
AFmatrix=data.frame()
for (i in levels(AFdata[,1]))
        {datanow<-AFdata[AFdata[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,6])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix=rbind(AFmatrix,wmatrix)
        }
names(AFmatrix)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Klet","AF_Kowa")


test<-as.character(AFdata$Chrom)
AFdata[,1]<-test
test<-as.character(AFmatrix$Chrom)
AFmatrix[,1]<-test

test2<-AFmatrix
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix[,i]))
	}
AFmatrix<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest", j, sep ="")
 	AF1<-AFmatrix[AFmatrix$Chrom==Candidates$Chr[j],]
	AFfreqofinterest<-AF1[AF1$Window_start>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_start<=(genemiddle[[1]][j]+windowsize[j]+100000)|AF1$Window_end>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_end<=(genemiddle[[1]][j]+windowsize[j]+100000),]
	assign(nam, AFfreqofinterest)
	}
rm(AFfreqofinterest)

SNPeffAFDabs<-SNPeffFst

AFdatasnp<-cbind(AFdata[,1:3],AFdata[,6])
names(AFdatasnp)<-c("CHROM","POS","AC","AC.1")
SNPeffdfAF<-merge(SNPeffAFDabs,AFdatasnp)
SNPeffdfAFo<-SNPeffdfAF[SNPeffdfAF$Effect=="MODERATE",]
SNPeffdfAFp<-SNPeffdfAF[SNPeffdfAF$Effect=="MODIFIER",]
SNPeffdfAFr<-SNPeffdfAF[SNPeffdfAF$Effect=="HIGH",]
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
	nam <- paste("AFDabs01percent_AFKletKowa_absoluteAFplots_candidate_genes_",j,"_",Candidates$Lyr_Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=18, height=18, units="cm", res=1000)
	plot(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS,ylab="AF difference Klet-Kowa",xlab="Scaffold position (bp)",ylim=c(0,1.4),cex.lab=1.5,xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]))
	rect(Candidates$gene_start[j],-2, Candidates$gene_end[j], 2,col="grey",border = NA)	
	points(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS)
	text(x=genemiddle[[1]][j],y=1.2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=1.1,x1=Candidates$gene_end[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Klet[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="darkblue")
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Kowa[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="lightblue")
	legend("top",legend=c("25SNP-window AF Klet","25SNP-window AF Kowa"),pch=15,col=c("darkblue","lightblue"),text.font=2,bg="white")
	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])],col="purple")
	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])],col="orange")
	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])],col="red")
	legend("topleft",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2)
	dev.off()
	}




############################################################################
#Candidate genes all metrics 0.1% lyrata genes without thaliana orthologues#
############################################################################
#combine genes from files in excel
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

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Chrom==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$POS>=(genemiddle[[1]][j]-windowsize[j])&AF1$POS<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")


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
rm(Genesofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

test<-as.factor(AFdata$Chrom)
AFdata[,1]<-test

winsize=25
AFmatrix=data.frame()
for (i in levels(AFdata[,1]))
        {datanow<-AFdata[AFdata[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,6])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix=rbind(AFmatrix,wmatrix)
        }
names(AFmatrix)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Klet","AF_Kowa")


test<-as.character(AFdata$Chrom)
AFdata[,1]<-test
test<-as.character(AFmatrix$Chrom)
AFmatrix[,1]<-test

test2<-AFmatrix
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix[,i]))
	}
AFmatrix<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest", j, sep ="")
 	AF1<-AFmatrix[AFmatrix$Chrom==Candidates$Chr[j],]
	AFfreqofinterest<-AF1[AF1$Window_start>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_start<=(genemiddle[[1]][j]+windowsize[j]+100000)|AF1$Window_end>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_end<=(genemiddle[[1]][j]+windowsize[j]+100000),]
	assign(nam, AFfreqofinterest)
	}
rm(AFfreqofinterest)

AFdatasnp<-cbind(AFdata[,1:3],AFdata[,6])
names(AFdatasnp)<-c("CHROM","POS","AC","AC.1")
SNPeffdfAF<-merge(SNPeffFst,AFdatasnp)
SNPeffdfAFo<-SNPeffdfAF[SNPeffdfAF$Effect=="MODERATE",]
SNPeffdfAFp<-SNPeffdfAF[SNPeffdfAF$Effect=="MODIFIER",]
SNPeffdfAFr<-SNPeffdfAF[SNPeffdfAF$Effect=="HIGH",]
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
	nam <- paste("Lyrataonly01percent_AFKletKowa_absoluteAFplots_candidate_genes_arenosa_",j,"_",Candidates$Lyr_Gene[j], sep ="")
	jpeg(paste(nam, '.jpeg', sep = ''), width=18, height=18, units="cm", res=1000)
	plot(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS,ylab="AF difference Klet-Kowa",xlab="Scaffold position (bp)",ylim=c(0,1.4),cex.lab=1.5,xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]))
	rect(Candidates$gene_start[j],-2, Candidates$gene_end[j], 2,col="grey",border = NA)	
	points(abs(get(paste("AFofinterest",j,sep=""))$AC-get(paste("AFofinterest",j,sep=""))$AC.1)~get(paste("AFofinterest",j,sep=""))$POS)
	text(x=genemiddle[[1]][j],y=1.2,labels=Candidates$Lyr_Gene[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=1.1,x1=Candidates$gene_end[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir2[i],length=0,col="red",lwd=3)
		}
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Klet[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="darkblue")
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Kowa[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="lightblue")
	legend("top",legend=c("25SNP-window AF Klet","25SNP-window AF Kowa"),pch=15,col=c("darkblue","lightblue"),text.font=2,bg="white")
	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Lyr_Gene[j])],col="purple")
	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Lyr_Gene[j])],col="orange")
	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Lyr_Gene[j])],col="red")
	legend("topleft",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2)
	dev.off()
	}



