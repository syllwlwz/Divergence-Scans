###########################
#Fay and Wu's H for lyrata#
###########################
data<-read.table("LyrataKletKowaGS_woGQfil_arenosa.csv",sep="\t",header=TRUE)
names(data)<-c("CHROM","POS","AC","CHROM.1","POS.1","AC.1")
data<-data.frame(data$CHROM,data$POS,data$AC,data$AC.1)

n1=36 #Klet
n2=32 #Kowa
nt=n1+n2
winsize=25

d<-data
names(d)<-c("CHROM","POS","AC","AC.1")

#whole data set

#Fay and Wu's H calculation

#exclude fixed sites within populations
d_Klet_1<-d[d$AC<max(d$AC),]
d_Klet<-d_Klet_1[,1:3]

d_Kowa_1<-d[d$AC.1<max(d$AC.1),]
d_Kowa<-data.frame(d_Kowa_1[,1:2],d_Kowa_1[,4])

#calculation of Pi
Pi_Klet=2*d_Klet[,3]*(n1-d_Klet[,3])/(n1*(n1-1))
Pi_Kowa=2*d_Kowa[,3]*(n2-d_Kowa[,3])/(n2*(n2-1))
Pi_Klet[1:3]
Pi_Kowa[1:3]

Theta_h_Klet=2*d_Klet[,3]*d_Klet[,3]/(n1*(n1-1))
Theta_h_Kowa=2*d_Kowa[,3]*d_Kowa[,3]/(n2*(n2-1))

#windows
#Klet
d_Klet_D<-data.frame(d_Klet,Pi_Klet,Theta_h_Klet)

FayandWusHdata=data.frame()
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
		Klet_Pi_win=mean(twin[,4])
		wmatrix[i,3]=Klet_Pi_win
		Klet_thetaH_win=mean(twin[,5])
		wmatrix[i,4]=Klet_thetaH_win
		wdt=wend+1
		wend=wend+winsize
		}
	wmatrix=cbind(scaff,wmatrix)
	FayandWusHdata<-rbind(FayandWusHdata,wmatrix)
	}
D=as.numeric(as.character(na.exclude(FayandWusHdata[,4])))-as.numeric(as.character(na.exclude(FayandWusHdata[,5])))
as.numeric(as.character(na.exclude(FayandWusHdata[,5])))[1:3]
D[1:3]
sd(D)
FayandWusH=D/sd(D)
FayandWusHdata<-na.exclude(FayandWusHdata)
FayandWusHdata[,4]<-FayandWusH
FayandWusHdata<-FayandWusHdata[,1:4]
write.table(FayandWusHdata,"FayandWusHKletlyrata_mean_arenosa.csv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

#Kowa
d_Kowa_D<-data.frame(d_Kowa,Pi_Kowa,Theta_h_Kowa)

FayandWusHdata=data.frame()
for (i in levels(d_Kowa_D[,1]))
        {datanow<-d_Kowa_D[d_Kowa_D[,1]==i,]
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
                Kowa_Pi_win=mean(twin[,4])
                wmatrix[i,3]=Kowa_Pi_win
                Kowa_thetaH_win=mean(twin[,5])
                wmatrix[i,4]=Kowa_thetaH_win
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        FayandWusHdata<-rbind(FayandWusHdata,wmatrix)
        }
D=as.numeric(as.character(na.exclude(FayandWusHdata[,4])))-as.numeric(as.character(na.exclude(FayandWusHdata[,5])))
D[1:3]
sd(D)
FayandWusH=D/sd(D)
FayandWusHdata<-na.exclude(FayandWusHdata)
FayandWusHdata[,4]<-FayandWusH
FayandWusHdata<-FayandWusHdata[,1:4]
write.table(FayandWusHdata,"FayandWusHKowalyrata_mean_arenosa.csv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


