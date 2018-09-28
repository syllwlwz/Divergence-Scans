##########################
#Nielsen 2dSFS for lyrata#
##########################
data<-read.table("LyrataKletKowaGS_woGQfil_arenosa.csv",sep="\t",header=TRUE)
names(data)<-c("CHROM","POS","AC","CHROM.1","POS.1","AC.1")

pos<-2
f1<-3
f2<-6
##your window size!
wsize=25
d<-data
names(d)<-c("CHROM","POS","AC","CHROM.1","POS.1","AC.1")
d=d[d$AC<max(d$AC) | d$AC.1<max(d$AC.1),]

d$AC[is.na(d[,3])]=0
d$AC.1[is.na(d[,6])]=0

Nielsendata<-data.frame()
count=0
#windows
for (i in levels(d[,1]))
	{
	count=count+1
	datanow<-d[d[,1]==i,]
	##initialize your matrix for 2 dimensional sfs
	m=matrix(nrow=length(min(datanow[,f1]):max(datanow[,f1])),ncol=length(min(datanow[,f2]):max(datanow[,f2])),data=0)
	##enter values for 2dSFS - rows are pop 1 sfs,  columns are pop 2 sfs
	for(i in min(datanow[,f1]):max(datanow[,f1])) 
		{
		if (min(datanow[,f1])=="0")
			{
			k=i+1
			for(j in min(datanow[,f2]):max(datanow[,f2])) 
				{
				if (min(datanow[,f2])=="0")
					{l=j+1
					m[k,l]=m[k,l]+nrow(datanow[datanow[,f1]==i & datanow[,f2]==j,])
					}
				else
					{l=which(min(datanow[,f2]):max(datanow[,f2])==j)
					m[k,l]=m[k,l]+nrow(datanow[datanow[,f1]==i & datanow[,f2]==j,])
					}
				}
			}
		else
			{
			k=which(min(datanow[,f1]):max(datanow[,f1])==i)
			for(j in min(datanow[,f2]):max(datanow[,f2])) 
				{
				if (min(datanow[,f2])=="0")
					{l=j+1
					m[k,l]=m[k,l]+nrow(datanow[datanow[,f1]==i & datanow[,f2]==j,])
					}
				else
					{l=which(min(datanow[,f2]):max(datanow[,f2])==j)
					m[k,l]=m[k,l]+nrow(datanow[datanow[,f1]==i & datanow[,f2]==j,])
					}
				}
			}

		}
	m[nrow(m),ncol(m)]=0 ##presumably you wont consider sites fixed in both pops?
	mn=m
	## convert number of sites to proportions (= genome-wide probabilities)
	m=mn/sum(mn)

	##calculate number of windows in your dataframe
	wcl=floor(length(datanow[,1])/wsize)
	s=1
	e=wsize
	##initialize matrix to record window info
	stat=matrix(0,wcl,5)
	for(w in 1:nrow(stat)) 
		{
		tw=datanow[s:e,]
		mw=matrix(nrow=length(min(datanow[,f1]):max(datanow[,f1])),ncol=length(min(datanow[,f2]):max(datanow[,f2])),data=0) ##Window-Specific Matrix!
		stat[w,1]=1
		stat[w,2]=mean(tw[,pos])
		stat[w,3]=min(tw[,pos])
		stat[w,4]=max(tw[,pos])
		LnL0=0
		LnL1=0
		if (min(datanow[,f1])=="0")
			{for(i in 0:max(datanow[,f1]))  ##Get 2dSFS for this window
				{
				k=i+1
				if (min(datanow[,f2])=="0")
					{for(j in 0:max(datanow[,f2])) 
						{
						l=j+1
						mw[k,l]=length(subset(tw[,f1],tw[,f1]==i & tw[,f2]==j))
						}
					}
				else
					{for(j in 1:max(datanow[,f2])) 
						{
						l=which(min(datanow[,f2]):max(datanow[,f2])==j)
						mw[k,l]=length(subset(tw[,f1],tw[,f1]==i & tw[,f2]==j))
						}
					}
				}
			}
		else
			{
			for(i in 1:max(datanow[,f1]))
				{
				k=which(min(datanow[,f1]):max(datanow[,f1])==i)
				if (min(datanow[,f2])=="0")
					{for(j in 0:max(datanow[,f2])) 
						{
						l=j+1
						mw[k,l]=length(subset(tw[,f1],tw[,f1]==i & tw[,f2]==j))
						}
					}
				else
					{for(j in 1:max(datanow[,f2])) 
						{
						l=which(min(datanow[,f2]):max(datanow[,f2])==j)
						mw[k,l]=length(subset(tw[,f1],tw[,f1]==i & tw[,f2]==j))
						}
					}
				}
			}
		mw[nrow(mw),ncol(mw)]=0
		pw=mw/sum(mw)
		LnL1=log(prod(pw^mw))
		if(LnL1==(-Inf)) 
			{
			LnL1=3^-600
			}
		LnL0=log(prod(m^mw))
		if(LnL0==(-Inf)) 
			{
			LnL0=3^-600
			}
		stat[w,5]=2*(LnL1-LnL0) ##your G-test statistic
		print(stat[w,5])
		s=e+1
		e=e+wsize
		}
	stat=cbind(count,stat)
	Nielsendata<-rbind(Nielsendata,stat)
	}

write.table(Nielsendata,"NielsenKletKowalyrata_arenosa.csv",quote=F,row.names=F,col.names=F)

