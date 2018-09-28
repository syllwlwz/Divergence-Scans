#########################
#FLK into 25 SNP windows#
#########################
data<-read.table("KletKowaFlk_alt.flk",sep=" ",header=TRUE)
winsize=25
d<-data[,2:6]

FLKdata<-data.frame()
scaff<-data.frame()
#windows
for (i in levels(d[,1]))
	{datanow<-d[d[,1]==i,]
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
		wmatrix[i,3]=min(twin[,5])
		wmatrix[i,4]=mean(twin[,4])
		wdt=wend+1
		wend=wend+winsize
		}
	wmatrix=cbind(scaff,wmatrix)
	FLKdata<-rbind(FLKdata,wmatrix)
	}
names(FLKdata)<-c("Scaffold","Window_start","Window_end","P-value","Flk")
FLKdata<-na.exclude(FLKdata)
write.table(FLKdata,"FlkKletKowalyrata_min_arenosa_alt.csv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
