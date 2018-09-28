data<-read.table("KletKowa_VarLD_complete.table",sep="\t",header=T)
w<-data.frame()
for (i in (levels(data[,1])))
	{d<-data[data[,1]==i,]
	datanew<-d[seq(1,nrow(d),25),]
	w<-rbind(w,datanew)
	}
levels(data[,1])
write.csv(w,"KletKowa_VarLD_25SNP_windows_arenosa",row.names=F,quote=F)
