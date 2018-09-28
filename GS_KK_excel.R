require(xlsx)


FSTmetlist<-read.csv("FST_KletKowa_01percent_arenosa.csv",header=T)
DXYmetlist<-read.csv("DXY_KletKowa_01percent_arenosa.csv",header=T)
AFDabsmetlist<-read.csv("AFDabs_KletKowa_01percent_arenosa.csv",header=T)
NIELSENmetlist<-read.csv("NIELSEN_KletKowa_01percent_arenosa.csv",header=T)
VARLDmetlist<-read.csv("VARLD_KletKowa_01percent_arenosa.csv",header=T)
FLKmetlist<-read.csv("FLK_KletKowa_01percent_arenosa.csv",header=T)
DDmetlist<-read.csv("DD_KletKowa_01percent_arenosa.csv",header=T)
HighSNPeffdiff<-read.table("GenesHIGHeffect_above40percentAFdiff_withUG_KK_arenosa.txt", sep="\t", header=T)


options(java.parameters = "-Xmx50000m")

write.xlsx2(FSTmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="Fst",col.names=TRUE,row.names=FALSE)
write.xlsx2(DXYmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="Dxy",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(AFDabsmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="AFDabs",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(NIELSENmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="Nielsen",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(VARLDmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="VarLD",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(FLKmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="Flk",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(DDmetlist,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="DD",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(HighSNPeffdiff,"Genes_01percent_KletKowa_reflyrata_arenosa.xlsx",sheetName="Highimpact",col.names=TRUE,row.names=FALSE,append=T)

