##################
#Post HC analysis#
##################

##################################################################
#Combining tables for both cohorts if not done by "paste" in unix#
##################################################################
data1<-read.table("Kowa_woGQfil_arenosa.table",header=TRUE,sep="\t")
data2<-read.table("Klet_woGQfil_arenosa.table",header=TRUE,sep="\t")
data<-data.frame(data2,data1)
str(data)
###################
#remove DP columns#
###################
LyrataKowaKletdata<-data.frame(data$CHROM,data$POS,data$AC,data$CHROM.1,data$POS.1,data$AC.1)
str(LyrataKowaKletdata)

#write.table(data,file="KowaKletLyrata_woGQfil.csv",row.name=FALSE)
#data<-read.table("KowaKletLyrata_woGQfil.csv",header=TRUE)
#LyrataKowaKletdata<-data[-c(4,5,9,10) ]

names(LyrataKowaKletdata)<-c("Chrom","POS","AC","Chrom","POS","AC")

write.table(LyrataKowaKletdata,"LyrataKletKowaGS_woGQfil_arenosa.csv",sep="\t",col.names=TRUE,row.names=FALSE)

