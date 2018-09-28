data<-read.table("KletKowa_GT.table",sep="\t",header=TRUE)

Klet_001_a07<-strsplit(as.character(data[,3]),"/")
Klet_001_a07_bind<-do.call(rbind,Klet_001_a07)
Klet_001_a08<-strsplit(as.character(data[,4]),"/")
Klet_001_a08_bind<-do.call(rbind,Klet_001_a08)
Klet_001_a09<-strsplit(as.character(data[,5]),"/")
Klet_001_a09_bind<-do.call(rbind,Klet_001_a09)
Klet_001_a10<-strsplit(as.character(data[,6]),"/")
Klet_001_a10_bind<-do.call(rbind,Klet_001_a10)
Klet_005_a18<-strsplit(as.character(data[,7]),"/")
Klet_005_a18_bind<-do.call(rbind,Klet_005_a18)
Klet_005_a19<-strsplit(as.character(data[,8]),"/")
Klet_005_a19_bind<-do.call(rbind,Klet_005_a19)
Klet_005_a21<-strsplit(as.character(data[,9]),"/")
Klet_005_a21_bind<-do.call(rbind,Klet_005_a21)
Klet_005_a28<-strsplit(as.character(data[,10]),"/")
Klet_005_a28_bind<-do.call(rbind,Klet_005_a28)
Klet_005_a32<-strsplit(as.character(data[,11]),"/")
Klet_005_a32_bind<-do.call(rbind,Klet_005_a32)
Kowa_001_a04<-strsplit(as.character(data[,12]),"/")
Kowa_001_a04_bind<-do.call(rbind,Kowa_001_a04)
Kowa_001_a05<-strsplit(as.character(data[,13]),"/")
Kowa_001_a05_bind<-do.call(rbind,Kowa_001_a05)
Kowa_001_a06<-strsplit(as.character(data[,14]),"/")
Kowa_001_a06_bind<-do.call(rbind,Kowa_001_a06)
Kowa_001_a07<-strsplit(as.character(data[,15]),"/")
Kowa_001_a07_bind<-do.call(rbind,Kowa_001_a07)
Kowa_001_a08<-strsplit(as.character(data[,16]),"/")
Kowa_001_a08_bind<-do.call(rbind,Kowa_001_a08)
Kowa_001_a09<-strsplit(as.character(data[,17]),"/")
Kowa_001_a09_bind<-do.call(rbind,Kowa_001_a09)
Kowa_001_a11<-strsplit(as.character(data[,18]),"/")
Kowa_001_a11_bind<-do.call(rbind,Kowa_001_a11)
Kowa_001_a12<-strsplit(as.character(data[,19]),"/")
Kowa_001_a12_bind<-do.call(rbind,Kowa_001_a12)

dataset<-cbind(Klet_001_a07_bind,Klet_001_a08_bind,Klet_001_a09_bind,Klet_001_a10_bind,Klet_005_a18_bind,Klet_005_a19_bind,Klet_005_a21_bind,Klet_005_a28_bind,Klet_005_a32_bind,Kowa_001_a04_bind,Kowa_001_a05_bind,Kowa_001_a06_bind,Kowa_001_a07_bind,Kowa_001_a08_bind,Kowa_001_a09_bind,Kowa_001_a11_bind,Kowa_001_a12_bind)
names(dataset)<-c("Klet_001_07_1","Klet_001_07_2","Klet_001_08_1","Klet_001_08_2","Klet_001_09_1","Klet_001_09_2","Klet_001_10_1","Klet_001_10_2","Klet_005_18_1","Klet_005_18_2","Klet_005_19_1","Klet_005_19_2","Klet_005_21_1","Klet_005_21_2","Klet_005_28_1","Klet_005_28_2","Klet_005_32_1","Klet_005_32_2","Kowa_001_04_1","Kowa_001_04_2","Kowa_001_05_1","Kowa_001_05_2","Kowa_001_06_1","Kowa_001_06_2","Kowa_001_07_1","Kowa_001_07_2","Kowa_001_08_1","Kowa_001_08_2","Kowa_001_09_1","Kowa_001_09_2","Kowa_001_11_1","Kowa_001_02_2","Kowa_001_12_1","Kowa_001_12_2")

Maxfrequ<-apply(dataset,1,max)

dataalt<-ifelse(dataset==Maxfrequ,"A","a")
#str(dataalt)

VarLDdata<-data.frame(data[,1:2])
for (i in seq(1,33,2))
	{datamerge1<-paste0(dataalt[,i],dataalt[,i+1])
	VarLDdata<-cbind(VarLDdata,datamerge1)
	}

VarLD1<-VarLDdata
VarLD1[]<-lapply(VarLD1,as.character)
VarLD1[VarLD1=="AA"]<-"1"
VarLD1[VarLD1=="Aa"]<-"2"
VarLD1[VarLD1=="aA"]<-"2"
VarLD1[VarLD1=="aa"]<-"3"
SNPID<-paste0(VarLD1[,1],"_",VarLD1[,2])
VarLD<-cbind(VarLD1[,1],SNPID,VarLD1[,2:19])
names(VarLD)<-NULL
write.table(VarLD,"VarLD_KletKowa_arenosa.csv",sep="\t",row.names=FALSE,col.names=FALSE,quote=F)
