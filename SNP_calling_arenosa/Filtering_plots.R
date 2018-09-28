data<-read.table("KletKowa_Filtest.table",header=TRUE,sep="\t")

library(ggplot2)

jpeg("QD_density.jpeg",width=26,height=18,units="cm",res=1000)
d <- density(na.exclude(data$QD))
plot(d)
dev.off()

jpeg("FS_density.jpeg",width=26,height=18,units="cm",res=1000)
FS<-density(na.exclude(data$FS))
plot(FS)
dev.off()

jpeg("FS_density_zoomin.jpeg",width=26,height=18,units="cm",res=1000)
plot(FS,xlim=c(0,60))
dev.off()

jpeg("SOR_density.jpeg",width=26,height=18,units="cm",res=1000)
SOR<-density(na.exclude(data$SOR))
plot(SOR)
dev.off()

jpeg("MQ_density.jpeg",width=26,height=18,units="cm",res=1000)
MQ<-density(na.exclude(data$MQ))
plot(MQ)
dev.off()

jpeg("MQ_density_zoomin.jpeg",width=26,height=18,units="cm",res=1000)
MQ<-density(na.exclude(data$MQ))
plot(MQ,xlim=c(20,70))
dev.off()

jpeg("MQRankSum_density.jpeg",width=26,height=18,units="cm",res=1000)
MQRankSum<-density(na.exclude(data$MQRankSum))
plot(MQRankSum)
dev.off()

jpeg("MQRankSum_density_zoomin.jpeg",width=26,height=18,units="cm",res=1000)
MQRankSum<-density(na.exclude(data$MQRankSum))
plot(MQRankSum,xlim=c(-15,10))
dev.off()

jpeg("ReadPosRankSum_density.jpeg",width=26,height=18,units="cm",res=1000)
ReadPosRankSum<-density(na.exclude(data$ReadPosRankSum))
plot(ReadPosRankSum)
dev.off()

jpeg("ReadPosRankSum_density_zoomin.jpeg",width=26,height=18,units="cm",res=1000)
ReadPosRankSum<-density(na.exclude(data$ReadPosRankSum))
plot(ReadPosRankSum,xlim=c(-8,8))
dev.off()


#    The read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.
#    Uninformative reads are not used in these annotations.


