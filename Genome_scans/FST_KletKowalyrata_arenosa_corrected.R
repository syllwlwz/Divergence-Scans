################
#FST for lyrata#
################
data<-read.table("LyrataKletKowaGS_woGQfil_arenosa.csv",sep="\t",header=TRUE)
names(data)<-c("CHROM","POS","AC","CHROM.1","POS.1","AC.1")
data<-data.frame(data$CHROM,data$POS,data$AC,data$AC.1)

n1=36 #Klet
n2=32 #Kowa
nt=n1+n2
winsize=25

d<-data
names(d)<-c("CHROM","POS","AC","AC.1")
d=d[d$AC<max(d$AC) | d$AC.1<max(d$AC.1),]

#whole data set
d$atot=as.numeric(d[,3])+as.numeric(d[,4])
d$ptot=d$atot/(nt)
d=subset(d,d$ptot<1)
d$ht=4*d$ptot^3*(1-d$ptot)+6*d$ptot^2*(1-d$ptot)^2+4*d$ptot*(1-d$ptot)^3
d[,3]=as.numeric(d[,3])/n1
d[,4]=as.numeric(d[,4])/n2
d$h1=4*d[,3]^3*(1-d[,3])+6*d[,3]^2*(1-d[,3])^2+4*d[,3]*(1-d[,3])^3
d$h2=4*d[,4]^3*(1-d[,4])+6*d[,4]^2*(1-d[,4])^2+4*d[,4]*(1-d[,4])^3
d$h12=((d$h1*n1)+(d$h2*n2))/nt
FSTdata<-data.frame()
scaff<-data.frame()
#windows
for (i in levels(d[,1]))
	{datanow<-d[d[,1]==i,]
	nwins=ceiling(nrow(datanow)/winsize)
	wmatrix=matrix(nrow=nwins,ncol=3,data=0)
	wdt=1
	wend=winsize
	scaff<-i
	for(i in 1:nrow(wmatrix)) 
		{
		twin=datanow[wdt:wend,]
		wmatrix[i,1]=min(twin[,2])
		wmatrix[i,2]=max(twin[,2])
		wmatrix[i,3]=mean(abs(twin$ht-twin$h12)/twin$ht)
		wdt=wend+1
		wend=wend+winsize
		}
	wmatrix=cbind(scaff,wmatrix)
	FSTdata<-rbind(FSTdata,wmatrix)
	}
write.table(FSTdata,"FstKletKowalyrata_arenosa_corrected.csv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

