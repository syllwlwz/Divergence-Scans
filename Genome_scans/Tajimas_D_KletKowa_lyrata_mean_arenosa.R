#######################
#Tajima's D for lyrata#
#######################
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

#Tajima's D calculation

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

#calculation of harmonic mean
a1=0
for (j in 1:(n1-1))
	{a1_temp<-1/j
	a1<-rbind(a1,a1_temp)
	}
a_Klet=sum(a1)
a_Klet
a2=0
for (j in 1:(n2-1))
	{a2_temp<-1/j
	a2<-rbind(a2,a2_temp)
	}
a_Kowa=sum(a2)
a_Kowa

#windows
#Klet
d_Klet_D<-data.frame(d_Klet,Pi_Klet)

Klet_theta_w_win=winsize/a_Klet/winsize
Klet_theta_w_win
TajimasDdata=data.frame()
for (i in levels(d_Klet_D[,1]))
	{datanow<-d_Klet_D[d_Klet_D[,1]==i,]
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
		Klet_Pi_win=mean(twin[,4])
		wmatrix[i,3]=Klet_Pi_win
		wdt=wend+1
		wend=wend+winsize
		}
	wmatrix=cbind(scaff,wmatrix)
	TajimasDdata<-rbind(TajimasDdata,wmatrix)
	}
D=as.numeric(as.character(na.exclude(TajimasDdata[,4])))-as.numeric(as.character(na.exclude(Klet_theta_w_win)))
D[1:3]
sd(D)
Dtest<-as.numeric(as.character(na.exclude(TajimasDdata[,4])))-as.numeric(as.character(na.exclude(Klet_theta_w_win/25)))
summary(Dtest)
Tajimas_D=D/sd(D)
TajimasDdata<-na.exclude(TajimasDdata)
TajimasDdata[,4]<-Tajimas_D
TajimasDdata[1:3,]
summary(TajimasDdata)
write.table(TajimasDdata,"TajimasDKletlyratamean_arenosa.csv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

#Kowa
d_Kowa_D<-data.frame(d_Kowa,Pi_Kowa)

Kowa_theta_w_win=winsize/a_Kowa/winsize
TajimasDdata=data.frame()
for (i in levels(d_Kowa_D[,1]))
	{datanow<-d_Kowa_D[d_Kowa_D[,1]==i,]
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
		Kowa_Pi_win=mean(twin[,4])
                wmatrix[i,3]=Kowa_Pi_win
		wdt=wend+1
		wend=wend+winsize
		}
	wmatrix=cbind(scaff,wmatrix)
	TajimasDdata<-rbind(TajimasDdata,wmatrix)
	}
D=as.numeric(as.character(na.exclude(TajimasDdata[,4])))-as.numeric(as.character(na.exclude(Kowa_theta_w_win)))
Tajimas_D=D/sd(D)
TajimasDdata<-na.exclude(TajimasDdata)
TajimasDdata[,4]<-Tajimas_D
write.table(TajimasDdata,"TajimasDKowalyratamean_arenosa.csv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


