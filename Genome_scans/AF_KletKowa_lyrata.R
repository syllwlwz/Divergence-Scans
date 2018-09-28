###############
#AF for lyrata#
###############
data<-read.table("LyrataKletKowaGS_woGQfil_arenosa.csv",sep="\t",header=TRUE)
names(data)<-c("CHROM","POS","AC","CHROM.1","POS.1","AC.1")
data<-data.frame(data$CHROM,data$POS,data$AC,data$AC.1)

n1=36 #Klet
n2=32 #Kowa
nt=n1+n2
winsize=25

d<-data
names(d)<-c("CHROM","POS","AC","AC.1")

#exclude fixed sites within populations
d_Klet_1<-d[d$AC<max(d$AC)|d$AC.1<max(d$AC.1),]
d_Klet_D<-d_Klet_1[,1:4]

Pimatrix=data.frame()
for (i in levels(d_Klet_D[,1]))
	{datanow<-d_Klet_D[d_Klet_D[,1]==i,]
	nwins=ceiling(nrow(datanow)/winsize)
	wmatrix=matrix(nrow=nwins,ncol=4,data=0)
	wdt=1
	wend=winsize
	scaff<-i
	for(i in 1:nrow(wmatrix)) 
		{
		twin=datanow[wdt:wend,]
		wmatrix[i,1]=min(twin[,2])
		wmatrix[i,2]=max(twin[,2])
		Klet_win=mean(twin[,3])
		wmatrix[i,3]=Klet_win
		wmatrix[i,4]=mean(twin[,4])
		wdt=wend+1
		wend=wend+winsize
		}
	wmatrix=cbind(scaff,wmatrix)
	Pimatrix=rbind(Pimatrix,wmatrix)
	}
write.table(na.exclude(Pimatrix),"AFKletKowalyratawindows_arenosa.csv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
