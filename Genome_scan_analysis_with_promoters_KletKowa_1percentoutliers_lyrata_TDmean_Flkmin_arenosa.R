######################################################################
#Analysis with potential promoter region before gene of 2 kb included#
######################################################################

#################################################
#load metrics and assign headers where necessary#
#################################################

fst<-read.table("FstKletKowalyrata_arenosa_corrected.csv",header=FALSE,sep="\t")
names(fst)<-c("Scaffold","Start_pos","End_pos","Fst")
fst<-subset(fst,!is.na(fst$Start_pos))

DDdata<-read.table("DDKletKowalyrata_arenosa.csv",sep="\t",header=F)
names(DDdata)<-c("Scaffold","Start_pos","End_pos","Mean_pos","Mean_allele_frequency_cohort1","Mean_allele_frequency_cohort2","Mean_abs_diff","Mean_pi_cohort1","Mean_pi_cohort2","Mean_raw_diff","DD")
DD<-data.frame(DDdata[,1:3],DDdata$DD)
names(DD)<-c("Scaffold","Start_pos","End_pos","DD")

AFD<-data.frame(DDdata[,1:3],DDdata$Mean_raw_diff)
names(AFD)<-c("Scaffold","Start_pos","End_pos","AFD")

AFDabs<-data.frame(DDdata[,1:3],DDdata$Mean_abs_diff)
names(AFDabs)<-c("Scaffold","Start_pos","End_pos","AFDabs")

Nielsen<-read.table("NielsenKletKowalyrata_arenosa.csv",header=FALSE)
Nielsen<-Nielsen[-c(2,3)]
names(Nielsen)<-c("Scaffold","Start_pos","End_pos","Nielsen_diversity")

Dxy<-read.table("DxyKletKowaarenosa.csv",header=FALSE,sep="\t")
Dxy<-subset(Dxy,!is.na(Dxy[,2]))
names(Dxy)<-c("Scaffold","Start_pos","End_pos","Dxy")

Flk<-read.table("FlkKletKowalyrata_min_arenosa_alt.csv",header=FALSE,sep="\t")
names(Flk)<-c("Scaffold","Start_pos","End_pos","pvalue","Flk")

VarLD<-read.table("KletKowa_VarLD_25SNP_windows_arenosa",header=T,sep=",")
VarLD<-VarLD[-c(2,3)]
names(VarLD)<-c("Scaffold","Position","VarLD")

TajimasD_Klet<-read.table("TajimasDKletlyratamean_arenosa.csv",header=FALSE,sep="\t")
names(TajimasD_Klet)<-c("Scaffold","Start_pos","End_pos","TajimasDKlet")

TajimasD_Kowa<-read.table("TajimasDKowalyratamean_arenosa.csv",header=FALSE,sep="\t")
names(TajimasD_Kowa)<-c("Scaffold","Start_pos","End_pos","TajimasDKowa")

FayandWusH_Klet<-read.table("FayandWusHKletlyrata_mean_arenosa.csv",header=FALSE,sep="\t")
names(FayandWusH_Klet)<-c("Scaffold","Start_pos","End_pos","FayandWusHKlet")

FayandWusH_Kowa<-read.table("FayandWusHKowalyrata_mean_arenosa.csv",header=FALSE,sep="\t")
names(FayandWusH_Kowa)<-c("Scaffold","Start_pos","End_pos","FayandWusHKowa")

AF<-read.table("AFKletKowalyratawindows_arenosa.csv",header=F,sep="\t")
names(AF)<-c("Scaffold","Start_pos","End_pos","AFKlet","AFKowa")
AF$AFKlet<-AF$AFKlet/36
AF$AFKowa<-AF$AFKowa/32

####################################
#combine all metrics into one table#
####################################

wholedata<-data.frame(fst,DD$DD,Nielsen$Nielsen_diversity,Dxy$Dxy,Flk$Flk,Flk$pvalue,VarLD$VarLD,AFD$AFD,AFDabs$AFD,DDdata$Mean_pi_cohort1,DDdata$Mean_pi_cohort2,AF$AFKlet,AF$AFKowa)
names(wholedata)<-c("Scaffold","Start_pos","End_pos","Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Klet","Pi_Kowa","AFKlet","AFKowa")
#write.table(wholedata,"AllpostUGtestsprewikill.csv",col.names=TRUE,sep="\t",row.names=FALSE)

#################################
#kill windows bigger than 100 kb#
#################################
wholedata<-wholedata[wholedata$End_pos-wholedata$Start_pos<100000,]
write.table(wholedata,"AllpostUGtests.csv",col.names=TRUE,sep="\t",row.names=FALSE)
#30 windows excluded
#191,957 windows left

#############################################
#calculate average and median of window size#
#############################################
mean(wholedata$End_pos-wholedata$Start_pos)
#858.9448
median(wholedata$End_pos-wholedata$Start_pos)
#380

#quantiles
quantile(wholedata$End_pos-wholedata$Start_pos,probs=c(0.01,0.1,0.25,0.75,0.9,0.99))
#      1%      10%      25%      75%      90%      99% 
#  137.00   211.00   275.00   582.00  1167.00 10666.44

#########################################
#find overlaps between windows and genes#
#########################################

library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)

#load gene lists
Lyr_ortholist=read.delim("LyrataGeneOrthogroup.txt", header=TRUE)
Thal_ortholist=read.delim("ThalianaGeneOrthogroup.txt", header=TRUE)
Ahal_ortholist=read.delim("HalleriGeneOrthogroup.txt", header=TRUE)

# Load the file of Thaliana descriptions
Thal_description=read.delim("TAIR10_functional_descriptions.txt", header=TRUE)
colnames(Thal_description)[1]<-"Gene"
test1<-substr(Thal_description[,1],1,9)
Thal_description[,1]<-test1

## Load GFF files and add promotor region to genes

### Step 1: import genes from GFF files

lyr_txdm = makeTxDbFromGFF("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt", format = c("gff3"), organism = "Arabidopsis lyrata")
#Ahal_txdm = makeTxDbFromGFF("annotation.gff", format = c("gff3"), organism = "Arabidopsis lyrata")

#View the information imported:
lyr_txdm
#Ahal_txdm

#Define a function that allows the addition of promotor regions
extend <- function(x, upstream=0, downstream=0) #from https://support.bioconductor.org/p/78652/
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

# Use the function to extend promotor region 2kb
  # Upstream=2000 adds 2kb upstream from transcription start site
  # Downstream=0 adds 0bp downstream from the gene.
#HalGenes_plusPromot = extend(genes(Ahal_txdm),upstream=2000, downstream=0)
LyrGenes_plusPromot = extend(genes(lyr_txdm),upstream=2000, downstream=0)

###########################
#select 1% outlier windows#
###########################

Perc1Fst<-wholedata[wholedata$Fst>=quantile(wholedata$Fst,0.99)[[1]],]
Perc1DD<-wholedata[wholedata$DD<=quantile(wholedata$DD,0.01)[[1]],]
Perc1Nielsen<-wholedata[wholedata$Nielsen>=quantile(wholedata$Nielsen,0.99)[[1]],]
Perc1Dxy<-wholedata[wholedata$Dxy>=quantile(wholedata$Dxy,0.99)[[1]],]
Perc1Flk<-wholedata[(wholedata$Flk>=quantile(wholedata$Flk,0.99)[[1]])&wholedata$Flk_pvalue<0.05,]
Perc1VarLD<-wholedata[wholedata$VarLD>=quantile(wholedata$VarLD,0.99)[[1]],]
Perc1AFD<-wholedata[wholedata$AFD>=quantile(wholedata$AFD,0.99)[[1]],]
Perc1AFDabs<-wholedata[wholedata$AFDabs>=quantile(wholedata$AFDabs,0.99)[[1]],]

### Convert list of significant markers to a "Genomic Ranges" table. 
Perc1Fst_GRange<-GRanges(seqnames=tolower(Perc1Fst$Scaffold),ranges=IRanges(start=Perc1Fst$Start_pos,end=Perc1Fst$End_pos))
Perc1DD_GRange<-GRanges(seqnames=tolower(Perc1DD$Scaffold),ranges=IRanges(start=Perc1DD$Start_pos,end=Perc1DD$End_pos))
Perc1Nielsen_GRange<-GRanges(seqnames=tolower(Perc1Nielsen$Scaffold),ranges=IRanges(start=Perc1Nielsen$Start_pos,end=Perc1Nielsen$End_pos))
Perc1DD_GRange<-GRanges(seqnames=tolower(Perc1DD$Scaffold),ranges=IRanges(start=Perc1DD$Start_pos,end=Perc1DD$End_pos))
Perc1Dxy_GRange<-GRanges(seqnames=tolower(Perc1Dxy$Scaffold),ranges=IRanges(start=Perc1Dxy$Start_pos,end=Perc1Dxy$End_pos))
Perc1Flk_GRange<-GRanges(seqnames=tolower(Perc1Flk$Scaffold),ranges=IRanges(start=Perc1Flk$Start_pos,end=Perc1Flk$End_pos))
Perc1VarLD_GRange<-GRanges(seqnames=tolower(Perc1VarLD$Scaffold),ranges=IRanges(start=Perc1VarLD$Start_pos,end=Perc1VarLD$End_pos))
Perc1AFD_GRange<-GRanges(seqnames=tolower(Perc1AFD$Scaffold),ranges=IRanges(start=Perc1AFD$Start_pos,end=Perc1AFD$End_pos))
Perc1AFDabs_GRange<-GRanges(seqnames=tolower(Perc1AFDabs$Scaffold),ranges=IRanges(start=Perc1AFDabs$Start_pos,end=Perc1AFDabs$End_pos))
values(Perc1Fst_GRange)<-Perc1Fst[,4:16]
values(Perc1DD_GRange)<-Perc1DD[,4:16]
values(Perc1Nielsen_GRange)<-Perc1Nielsen[,4:16]
values(Perc1DD_GRange)<-Perc1DD[,4:16]
values(Perc1Dxy_GRange)<-Perc1Dxy[,4:16]
values(Perc1Flk_GRange)<-Perc1Flk[,4:16]
values(Perc1VarLD_GRange)<-Perc1VarLD[,4:16]
values(Perc1AFD_GRange)<-Perc1AFD[,4:16]
values(Perc1AFDabs_GRange)<-Perc1AFDabs[,4:16]

## Merge Significant Regions with Gene list in the lyrata Genome

Perc1Fst_lyratagenes=mergeByOverlaps(Perc1Fst_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1Fst_lyratagenes_df=data.frame(as.data.frame(Perc1Fst_lyratagenes$Perc1Fst_GRange),as.data.frame(Perc1Fst_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1Fst_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Klet","Pi_Kowa","AF_Klet","AF_Kowa","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1Fst_lyratagenes_df$gene_start<-ifelse(Perc1Fst_lyratagenes_df$gene_strand=="+",Perc1Fst_lyratagenes_df$gene_start+2000,Perc1Fst_lyratagenes_df$gene_start)
Perc1Fst_lyratagenes_df$gene_end<-ifelse(Perc1Fst_lyratagenes_df$gene_strand=="-",Perc1Fst_lyratagenes_df$gene_end-2000,Perc1Fst_lyratagenes_df$gene_end)
Perc1Fst_lyratagenes_df$gene_size<-Perc1Fst_lyratagenes_df$gene_size-2000

Perc1DD_lyratagenes=mergeByOverlaps(Perc1DD_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1DD_lyratagenes_df=data.frame(as.data.frame(Perc1DD_lyratagenes$Perc1DD_GRange),as.data.frame(Perc1DD_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1DD_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Klet","Pi_Kowa","AF_Klet","AF_Kowa","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1DD_lyratagenes_df$gene_start<-ifelse(Perc1DD_lyratagenes_df$gene_strand=="+",Perc1DD_lyratagenes_df$gene_start+2000,Perc1DD_lyratagenes_df$gene_start)
Perc1DD_lyratagenes_df$gene_end<-ifelse(Perc1DD_lyratagenes_df$gene_strand=="-",Perc1DD_lyratagenes_df$gene_end-2000,Perc1DD_lyratagenes_df$gene_end)
Perc1DD_lyratagenes_df$gene_size<-Perc1DD_lyratagenes_df$gene_size-2000

Perc1Nielsen_lyratagenes=mergeByOverlaps(Perc1Nielsen_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1Nielsen_lyratagenes_df=data.frame(as.data.frame(Perc1Nielsen_lyratagenes$Perc1Nielsen_GRange),as.data.frame(Perc1Nielsen_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1Nielsen_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Klet","Pi_Kowa","AF_Klet","AF_Kowa","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1Nielsen_lyratagenes_df$gene_start<-ifelse(Perc1Nielsen_lyratagenes_df$gene_strand=="+",Perc1Nielsen_lyratagenes_df$gene_start+2000,Perc1Nielsen_lyratagenes_df$gene_start)
Perc1Nielsen_lyratagenes_df$gene_end<-ifelse(Perc1Nielsen_lyratagenes_df$gene_strand=="-",Perc1Nielsen_lyratagenes_df$gene_end-2000,Perc1Nielsen_lyratagenes_df$gene_end)
Perc1Nielsen_lyratagenes_df$gene_size<-Perc1Nielsen_lyratagenes_df$gene_size-2000

Perc1Dxy_lyratagenes=mergeByOverlaps(Perc1Dxy_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1Dxy_lyratagenes_df=data.frame(as.data.frame(Perc1Dxy_lyratagenes$Perc1Dxy_GRange),as.data.frame(Perc1Dxy_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1Dxy_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Klet","Pi_Kowa","AF_Klet","AF_Kowa", "Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1Dxy_lyratagenes_df$gene_start<-ifelse(Perc1Dxy_lyratagenes_df$gene_strand=="+",Perc1Dxy_lyratagenes_df$gene_start+2000,Perc1Dxy_lyratagenes_df$gene_start)
Perc1Dxy_lyratagenes_df$gene_end<-ifelse(Perc1Dxy_lyratagenes_df$gene_strand=="-",Perc1Dxy_lyratagenes_df$gene_end-2000,Perc1Dxy_lyratagenes_df$gene_end)
Perc1Dxy_lyratagenes_df$gene_size<-Perc1Dxy_lyratagenes_df$gene_size-2000

Perc1Flk_lyratagenes=mergeByOverlaps(Perc1Flk_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1Flk_lyratagenes_df=data.frame(as.data.frame(Perc1Flk_lyratagenes$Perc1Flk_GRange),as.data.frame(Perc1Flk_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1Flk_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Klet","Pi_Kowa","AF_Klet","AF_Kowa","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1Flk_lyratagenes_df$gene_start<-ifelse(Perc1Flk_lyratagenes_df$gene_strand=="+",Perc1Flk_lyratagenes_df$gene_start+2000,Perc1Flk_lyratagenes_df$gene_start)
Perc1Flk_lyratagenes_df$gene_end<-ifelse(Perc1Flk_lyratagenes_df$gene_strand=="-",Perc1Flk_lyratagenes_df$gene_end-2000,Perc1Flk_lyratagenes_df$gene_end)
Perc1Flk_lyratagenes_df$gene_size<-Perc1Flk_lyratagenes_df$gene_size-2000

Perc1VarLD_lyratagenes=mergeByOverlaps(Perc1VarLD_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1VarLD_lyratagenes_df=data.frame(as.data.frame(Perc1VarLD_lyratagenes$Perc1VarLD_GRange),as.data.frame(Perc1VarLD_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1VarLD_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Klet","Pi_Kowa","AF_Klet","AF_Kowa","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1VarLD_lyratagenes_df$gene_start<-ifelse(Perc1VarLD_lyratagenes_df$gene_strand=="+",Perc1VarLD_lyratagenes_df$gene_start+2000,Perc1VarLD_lyratagenes_df$gene_start)
Perc1VarLD_lyratagenes_df$gene_end<-ifelse(Perc1VarLD_lyratagenes_df$gene_strand=="-",Perc1VarLD_lyratagenes_df$gene_end-2000,Perc1VarLD_lyratagenes_df$gene_end)
Perc1VarLD_lyratagenes_df$gene_size<-Perc1VarLD_lyratagenes_df$gene_size-2000

Perc1AFD_lyratagenes=mergeByOverlaps(Perc1AFD_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1AFD_lyratagenes_df=data.frame(as.data.frame(Perc1AFD_lyratagenes$Perc1AFD_GRange),as.data.frame(Perc1AFD_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1AFD_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Klet","Pi_Kowa","AF_Klet","AF_Kowa","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1AFD_lyratagenes_df$gene_start<-ifelse(Perc1AFD_lyratagenes_df$gene_strand=="+",Perc1AFD_lyratagenes_df$gene_start+2000,Perc1AFD_lyratagenes_df$gene_start)
Perc1AFD_lyratagenes_df$gene_end<-ifelse(Perc1AFD_lyratagenes_df$gene_strand=="-",Perc1AFD_lyratagenes_df$gene_end-2000,Perc1AFD_lyratagenes_df$gene_end)
Perc1AFD_lyratagenes_df$gene_size<-Perc1AFD_lyratagenes_df$gene_size-2000

Perc1AFDabs_lyratagenes=mergeByOverlaps(Perc1AFDabs_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1AFDabs_lyratagenes_df=data.frame(as.data.frame(Perc1AFDabs_lyratagenes$Perc1AFDabs_GRange),as.data.frame(Perc1AFDabs_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1AFDabs_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Klet","Pi_Kowa","AF_Klet","AF_Kowa","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1AFDabs_lyratagenes_df$gene_start<-ifelse(Perc1AFDabs_lyratagenes_df$gene_strand=="+",Perc1AFDabs_lyratagenes_df$gene_start+2000,Perc1AFDabs_lyratagenes_df$gene_start)
Perc1AFDabs_lyratagenes_df$gene_end<-ifelse(Perc1AFDabs_lyratagenes_df$gene_strand=="-",Perc1AFDabs_lyratagenes_df$gene_end-2000,Perc1AFDabs_lyratagenes_df$gene_end)
Perc1AFDabs_lyratagenes_df$gene_size<-Perc1AFDabs_lyratagenes_df$gene_size-2000


######                                                                     #####
###### REGIONS NOT OVERLAPPING GENES WILL BE LOST FROM LIST AT THIS POINT  #####
######                                                                     #####

### Create a file with the list of Thaliana genes and descriptions

#lyrata genes to Orthogroup
Perc1Fst_lyratagenes_df_w_orthogroup= merge(Perc1Fst_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1DD_lyratagenes_df_w_orthogroup= merge(Perc1DD_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1Nielsen_lyratagenes_df_w_orthogroup= merge(Perc1Nielsen_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1DD_lyratagenes_df_w_orthogroup= merge(Perc1DD_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1Dxy_lyratagenes_df_w_orthogroup= merge(Perc1Dxy_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1Flk_lyratagenes_df_w_orthogroup= merge(Perc1Flk_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1VarLD_lyratagenes_df_w_orthogroup= merge(Perc1VarLD_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1AFD_lyratagenes_df_w_orthogroup= merge(Perc1AFD_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1AFDabs_lyratagenes_df_w_orthogroup= merge(Perc1AFDabs_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)

#Orthogroup to Thaliana Genes
Perc1Fst_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1Fst_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1Fst_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1Fst_lyratagenes_df_w_orthogroup_Thalgenes)[26]="Gene"

Perc1DD_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1DD_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1DD_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1DD_lyratagenes_df_w_orthogroup_Thalgenes)[26]="Gene"

Perc1Nielsen_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1Nielsen_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1Nielsen_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1Nielsen_lyratagenes_df_w_orthogroup_Thalgenes)[26]="Gene"

Perc1Dxy_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1Dxy_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1Dxy_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1Dxy_lyratagenes_df_w_orthogroup_Thalgenes)[26]="Gene"

Perc1Flk_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1Flk_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1Flk_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1Flk_lyratagenes_df_w_orthogroup_Thalgenes)[26]="Gene"

Perc1VarLD_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1VarLD_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1VarLD_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1VarLD_lyratagenes_df_w_orthogroup_Thalgenes)[26]="Gene"

Perc1AFD_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1AFD_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1AFD_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1AFD_lyratagenes_df_w_orthogroup_Thalgenes)[26]="Gene"

Perc1AFDabs_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1AFDabs_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1AFDabs_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1AFDabs_lyratagenes_df_w_orthogroup_Thalgenes)[26]="Gene"


#Thaliana Genes to Gene Descriptions
Thalgenes_described_Perc1Fst_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1Fst_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Thalgenes_described_Perc1DD_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1DD_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Thalgenes_described_Perc1Nielsen_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1Nielsen_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Thalgenes_described_Perc1Dxy_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1Dxy_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Thalgenes_described_Perc1Flk_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1Flk_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Thalgenes_described_Perc1VarLD_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1VarLD_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Thalgenes_described_Perc1AFD_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1AFD_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Thalgenes_described_Perc1AFDabs_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1AFDabs_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)

#Generate file
write.table(Thalgenes_described_Perc1Fst_lyratagenes_df_w_orthogroup,"Genes1percentFst_lyratagenome_TDmean_arenosa.txt", sep="\t", row.names=F)
write.table(Thalgenes_described_Perc1DD_lyratagenes_df_w_orthogroup,"Genes1percentDD_lyratagenome_TDmean_arenosa.txt", sep="\t", row.names=F)
write.table(Thalgenes_described_Perc1Nielsen_lyratagenes_df_w_orthogroup,"Genes1percentNielsen_lyratagenome_TDmean_arenosa.txt", sep="\t", row.names=F)
write.table(Thalgenes_described_Perc1Dxy_lyratagenes_df_w_orthogroup,"Genes1percentDxy_lyratagenome_TDmean_arenosa.txt", sep="\t", row.names=F)
write.table(Thalgenes_described_Perc1Flk_lyratagenes_df_w_orthogroup,"Genes1percentFlk_lyratagenome_TDmean_arenosa.txt", sep="\t", row.names=F)
write.table(Thalgenes_described_Perc1VarLD_lyratagenes_df_w_orthogroup,"Genes1percentVarLD_lyratagenome_TDmean_arenosa.txt", sep="\t", row.names=F)
write.table(Thalgenes_described_Perc1AFD_lyratagenes_df_w_orthogroup,"Genes1percentAFD_lyratagenome_TDmean_arenosa.txt", sep="\t", row.names=F)
write.table(Thalgenes_described_Perc1AFDabs_lyratagenes_df_w_orthogroup,"Genes1percentAFD_lyratagenome_TDmean_arenosa.txt", sep="\t", row.names=F)

###############################
#Tajima's D and Fay and Wu's H#
###############################

####################################
#combine all metrics into one table#
####################################

TDandFWH_Klet<-data.frame(TajimasD_Klet,FayandWusH_Klet$FayandWusHKlet)
names(TDandFWH_Klet)<-c("Scaffold","Start_pos","End_pos","TajimasD_Klet","FayandWusH_Klet")

TDandFWH_Kowa<-data.frame(TajimasD_Kowa,FayandWusH_Kowa$FayandWusHKowa)
names(TDandFWH_Kowa)<-c("Scaffold","Start_pos","End_pos","TajimasD_Kowa","FayandWusH_Kowa")

#189,291 windows Klet, 187,178 windows Kowa
#write.table(wholedata,"AllpostUGtestsprewikill.csv",col.names=TRUE,sep="\t",row.names=FALSE)

#################################
#kill windows bigger than 100 kb#
#################################
TDandFWH_Klet<-TDandFWH_Klet[TDandFWH_Klet$End_pos-TDandFWH_Klet$Start_pos<100000,]
TDandFWH_Kowa<-TDandFWH_Kowa[TDandFWH_Kowa$End_pos-TDandFWH_Kowa$Start_pos<100000,]
write.table(TDandFWH_Klet,"TDFWH_Klet.csv",col.names=TRUE,sep="\t",row.names=FALSE)
write.table(TDandFWH_Kowa,"TDFWH_Kowa.csv",col.names=TRUE,sep="\t",row.names=FALSE)
#189,263 windows Klet, 187,151 windows Kowa)

#############################################
#calculate average and median of window size#
#############################################
mean(TDandFWH_Klet$End_pos-TDandFWH_Klet$Start_pos)
#870.8487
median(TDandFWH_Klet$End_pos-TDandFWH_Klet$Start_pos)
#386

#quantiles
quantile(TDandFWH_Klet$End_pos-TDandFWH_Klet$Start_pos,probs=c(0.01,0.1,0.25,0.75,0.9,0.99))
#      1%      10%      25%      75%      90%      99% 
#  138.00   213.00   278.00   594.00  1197.00 10838.66

mean(TDandFWH_Kowa$End_pos-TDandFWH_Kowa$Start_pos)
#883.4407
median(TDandFWH_Kowa$End_pos-TDandFWH_Kowa$Start_pos)
#390

#quantiles
quantile(TDandFWH_Kowa$End_pos-TDandFWH_Kowa$Start_pos,probs=c(0.01,0.1,0.25,0.75,0.9,0.99))
#     1%     10%     25%     75%     90%     99% 
#  138.0   215.0   281.0   603.0  1216.0 10980.5

Perc1TajimasD_Klet<-TDandFWH_Klet[TDandFWH_Klet$TajimasD_Klet<=quantile(TDandFWH_Klet$TajimasD_Klet,0.01)[[1]],]
Perc1TajimasD_Kowa<-TDandFWH_Kowa[TDandFWH_Kowa$TajimasD_Kowa<=quantile(TDandFWH_Kowa$TajimasD_Kowa,0.01)[[1]],]
Perc1FayandWusH_Klet<-TDandFWH_Klet[TDandFWH_Klet$FayandWusH_Klet<=quantile(TDandFWH_Klet$FayandWusH_Klet,0.01)[[1]],]
Perc1FayandWusH_Kowa<-TDandFWH_Kowa[TDandFWH_Kowa$FayandWusH_Kowa<=quantile(TDandFWH_Kowa$FayandWusH_Kowa,0.01)[[1]],]

Perc1TajimasD_Klet_GRange<-GRanges(seqnames=tolower(Perc1TajimasD_Klet$Scaffold),ranges=IRanges(start=Perc1TajimasD_Klet$Start_pos,end=Perc1TajimasD_Klet$End_pos))
Perc1TajimasD_Kowa_GRange<-GRanges(seqnames=tolower(Perc1TajimasD_Kowa$Scaffold),ranges=IRanges(start=Perc1TajimasD_Kowa$Start_pos,end=Perc1TajimasD_Kowa$End_pos))
Perc1FayandWusH_Klet_GRange<-GRanges(seqnames=tolower(Perc1FayandWusH_Klet$Scaffold),ranges=IRanges(start=Perc1FayandWusH_Klet$Start_pos,end=Perc1FayandWusH_Klet$End_pos))
Perc1FayandWusH_Kowa_GRange<-GRanges(seqnames=tolower(Perc1FayandWusH_Kowa$Scaffold),ranges=IRanges(start=Perc1FayandWusH_Kowa$Start_pos,end=Perc1FayandWusH_Kowa$End_pos))

values(Perc1TajimasD_Klet_GRange)<-Perc1TajimasD_Klet[,4:5]
values(Perc1TajimasD_Kowa_GRange)<-Perc1TajimasD_Kowa[,4:5]
values(Perc1FayandWusH_Klet_GRange)<-Perc1FayandWusH_Klet[,4:5]
values(Perc1FayandWusH_Kowa_GRange)<-Perc1FayandWusH_Kowa[,4:5]

## Merge Significant Regions with Gene list in the Halleri Genome

Perc1TajimasD_Klet_lyratagenes=mergeByOverlaps(Perc1TajimasD_Klet_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1TajimasD_Klet_lyratagenes_df=data.frame(as.data.frame(Perc1TajimasD_Klet_lyratagenes$Perc1TajimasD_Klet_GRange),as.data.frame(Perc1TajimasD_Klet_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1TajimasD_Klet_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "TajimasD_Klet","FayandWusH_Klet", "Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1TajimasD_Klet_lyratagenes_df$gene_start<-ifelse(Perc1TajimasD_Klet_lyratagenes_df$gene_strand=="+",Perc1TajimasD_Klet_lyratagenes_df$gene_start+2000,Perc1TajimasD_Klet_lyratagenes_df$gene_start)
Perc1TajimasD_Klet_lyratagenes_df$gene_end<-ifelse(Perc1TajimasD_Klet_lyratagenes_df$gene_strand=="-",Perc1TajimasD_Klet_lyratagenes_df$gene_end-2000,Perc1TajimasD_Klet_lyratagenes_df$gene_end)
Perc1TajimasD_Klet_lyratagenes_df$gene_size<-Perc1TajimasD_Klet_lyratagenes_df$gene_size-2000

Perc1TajimasD_Kowa_lyratagenes=mergeByOverlaps(Perc1TajimasD_Kowa_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1TajimasD_Kowa_lyratagenes_df=data.frame(as.data.frame(Perc1TajimasD_Kowa_lyratagenes$Perc1TajimasD_Kowa_GRange),as.data.frame(Perc1TajimasD_Kowa_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1TajimasD_Kowa_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "TajimasD_Kowa","FayandWusH_Kowa", "Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1TajimasD_Kowa_lyratagenes_df$gene_start<-ifelse(Perc1TajimasD_Kowa_lyratagenes_df$gene_strand=="+",Perc1TajimasD_Kowa_lyratagenes_df$gene_start+2000,Perc1TajimasD_Kowa_lyratagenes_df$gene_start)
Perc1TajimasD_Kowa_lyratagenes_df$gene_end<-ifelse(Perc1TajimasD_Kowa_lyratagenes_df$gene_strand=="-",Perc1TajimasD_Kowa_lyratagenes_df$gene_end-2000,Perc1TajimasD_Kowa_lyratagenes_df$gene_end)
Perc1TajimasD_Kowa_lyratagenes_df$gene_size<-Perc1TajimasD_Kowa_lyratagenes_df$gene_size-2000

Perc1FayandWusH_Klet_lyratagenes=mergeByOverlaps(Perc1FayandWusH_Klet_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1FayandWusH_Klet_lyratagenes_df=data.frame(as.data.frame(Perc1FayandWusH_Klet_lyratagenes$Perc1FayandWusH_Klet_GRange),as.data.frame(Perc1FayandWusH_Klet_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1FayandWusH_Klet_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "TajimasD_Klet","FayandWusH_Klet", "Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1FayandWusH_Klet_lyratagenes_df$gene_start<-ifelse(Perc1FayandWusH_Klet_lyratagenes_df$gene_strand=="+",Perc1FayandWusH_Klet_lyratagenes_df$gene_start+2000,Perc1FayandWusH_Klet_lyratagenes_df$gene_start)
Perc1FayandWusH_Klet_lyratagenes_df$gene_end<-ifelse(Perc1FayandWusH_Klet_lyratagenes_df$gene_strand=="-",Perc1FayandWusH_Klet_lyratagenes_df$gene_end-2000,Perc1FayandWusH_Klet_lyratagenes_df$gene_end)
Perc1FayandWusH_Klet_lyratagenes_df$gene_size<-Perc1FayandWusH_Klet_lyratagenes_df$gene_size-2000

Perc1FayandWusH_Kowa_lyratagenes=mergeByOverlaps(Perc1FayandWusH_Kowa_GRange,LyrGenes_plusPromot,type=c("any"))
Perc1FayandWusH_Kowa_lyratagenes_df=data.frame(as.data.frame(Perc1FayandWusH_Kowa_lyratagenes$Perc1FayandWusH_Kowa_GRange),as.data.frame(Perc1FayandWusH_Kowa_lyratagenes$LyrGenes_plusPromot))
colnames(Perc1FayandWusH_Kowa_lyratagenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "TajimasD_Kowa","FayandWusH_Kowa", "Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1FayandWusH_Kowa_lyratagenes_df$gene_start<-ifelse(Perc1FayandWusH_Kowa_lyratagenes_df$gene_strand=="+",Perc1FayandWusH_Kowa_lyratagenes_df$gene_start+2000,Perc1FayandWusH_Kowa_lyratagenes_df$gene_start)
Perc1FayandWusH_Kowa_lyratagenes_df$gene_end<-ifelse(Perc1FayandWusH_Kowa_lyratagenes_df$gene_strand=="-",Perc1FayandWusH_Kowa_lyratagenes_df$gene_end-2000,Perc1FayandWusH_Kowa_lyratagenes_df$gene_end)
Perc1FayandWusH_Kowa_lyratagenes_df$gene_size<-Perc1FayandWusH_Kowa_lyratagenes_df$gene_size-2000

######                                                                     #####
###### REGIONS NOT OVERLAPPING GENES WILL BE LOST FROM LIST AT THIS POINT  #####
######                                                                     #####

### Create a file with the list of Thaliana genes and descriptions

#lyrata genes to Orthogroup
Perc1TajimasD_Klet_lyratagenes_df_w_orthogroup= merge(Perc1TajimasD_Klet_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1TajimasD_Kowa_lyratagenes_df_w_orthogroup= merge(Perc1TajimasD_Kowa_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1FayandWusH_Klet_lyratagenes_df_w_orthogroup= merge(Perc1FayandWusH_Klet_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)
Perc1FayandWusH_Kowa_lyratagenes_df_w_orthogroup= merge(Perc1FayandWusH_Kowa_lyratagenes_df, Lyr_ortholist, by="Gene", all.x=TRUE)

#Orthogroup to Thaliana Genes
Perc1TajimasD_Klet_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1TajimasD_Klet_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1TajimasD_Klet_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1TajimasD_Klet_lyratagenes_df_w_orthogroup_Thalgenes)[15]="Gene"

Perc1TajimasD_Kowa_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1TajimasD_Kowa_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1TajimasD_Kowa_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1TajimasD_Kowa_lyratagenes_df_w_orthogroup_Thalgenes)[15]="Gene"

Perc1FayandWusH_Klet_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1FayandWusH_Klet_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1FayandWusH_Klet_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1FayandWusH_Klet_lyratagenes_df_w_orthogroup_Thalgenes)[15]="Gene"

Perc1FayandWusH_Kowa_lyratagenes_df_w_orthogroup_Thalgenes= merge(Perc1FayandWusH_Kowa_lyratagenes_df_w_orthogroup, Thal_ortholist, by="OrthoGroup", all.x=TRUE)
colnames(Perc1FayandWusH_Kowa_lyratagenes_df_w_orthogroup_Thalgenes)[2]="Lyr_Gene"
colnames(Perc1FayandWusH_Kowa_lyratagenes_df_w_orthogroup_Thalgenes)[15]="Gene"

#Thaliana Genes to Gene Descriptions
Thalgenes_described_Perc1TajimasD_Klet_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1TajimasD_Klet_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Thalgenes_described_Perc1TajimasD_Kowa_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1TajimasD_Kowa_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Thalgenes_described_Perc1FayandWusH_Klet_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1FayandWusH_Klet_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)
Thalgenes_described_Perc1FayandWusH_Kowa_lyratagenes_df_w_orthogroup= merge(Thal_description, Perc1FayandWusH_Kowa_lyratagenes_df_w_orthogroup_Thalgenes, by ="Gene",all.y=T)

#Generate file
write.table(Thalgenes_described_Perc1TajimasD_Klet_lyratagenes_df_w_orthogroup,"Genes1percentTajimasD_Klet_lyratagenome_TDmean.txt", sep="\t", row.names=F)
write.table(Thalgenes_described_Perc1TajimasD_Kowa_lyratagenes_df_w_orthogroup,"Genes1percentTajimasD_Kowa_lyratagenome_TDmean.txt", sep="\t", row.names=F)
write.table(Thalgenes_described_Perc1FayandWusH_Klet_lyratagenes_df_w_orthogroup,"Genes1percentFayandWusH_Klet_lyratagenome_TDmean.txt", sep="\t", row.names=F)
write.table(Thalgenes_described_Perc1FayandWusH_Kowa_lyratagenes_df_w_orthogroup,"Genes1percentFayandWusH_Kowa_lyratagenome_TDmean.txt", sep="\t", row.names=F)

options(java.parameters = "-Xmx12000m")

require(xlsx)
write.xlsx2(Thalgenes_described_Perc1Fst_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="Fst_1%",col.names=TRUE,row.names=FALSE)
write.xlsx2(Thalgenes_described_Perc1DD_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="DD_1%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Thalgenes_described_Perc1Nielsen_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="Nielsen_1%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Thalgenes_described_Perc1Dxy_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="Dxy_1%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Thalgenes_described_Perc1Flk_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="Flk_1%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Thalgenes_described_Perc1VarLD_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="VarLD_1%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Thalgenes_described_Perc1AFD_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="AFD_1%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Thalgenes_described_Perc1AFDabs_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="AFDabs_1%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Thalgenes_described_Perc1TajimasD_Klet_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="TajimasD_Klet_1%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Thalgenes_described_Perc1TajimasD_Kowa_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="TajimasD_Kowa_1%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Thalgenes_described_Perc1FayandWusH_Klet_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="FayandWusH_Klet_1%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Thalgenes_described_Perc1FayandWusH_Kowa_lyratagenes_df_w_orthogroup,"GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx",sheetName="FayandWusH_Kowa_1%",col.names=TRUE,row.names=FALSE,append=TRUE)


#select values for % outliers#
FSTPerc<-rbind(names(quantile(wholedata$Fst,probs = seq(0.99,1,0.001))),quantile(wholedata$Fst,probs = seq(0.99,1,0.001)))
AFDabsPerc<-rbind(names(quantile(wholedata$AFDabs,probs = seq(0.99,1,0.001))),quantile(wholedata$AFDabs,probs = seq(0.99,1,0.001)))
NIELSENPerc<-rbind(names(quantile(wholedata$Nielsen,probs = seq(0.99,1,0.001))),quantile(wholedata$Nielsen,probs = seq(0.99,1,0.001)))
DXYPerc<-rbind(names(quantile(wholedata$Dxy,probs = seq(0.99,1,0.001))),quantile(wholedata$Dxy,probs = seq(0.99,1,0.001)))
VARLDPerc<-rbind(names(quantile(wholedata$VarLD,probs = seq(0.99,1,0.001))),quantile(wholedata$VarLD,probs = seq(0.99,1,0.001)))
FLKPerc<-rbind(names(quantile(wholedata$Flk,probs = seq(0.99,1,0.001))),quantile(wholedata$Flk,probs = seq(0.99,1,0.001)))
TDPercKlet<-rbind(names(quantile(TDandFWH_Klet$TajimasD_Klet,probs = seq(0,0.01,0.001))),quantile(TDandFWH_Klet$TajimasD_Klet,probs = seq(0,0.01,0.001)))
TDPercKowa<-rbind(names(quantile(TDandFWH_Kowa$TajimasD_Kowa,probs = seq(0,0.01,0.001))),quantile(TDandFWH_Kowa$TajimasD_Kowa,probs = seq(0,0.01,0.001)))
FWHPercKlet<-rbind(names(quantile(TDandFWH_Klet$FayandWusH_Klet,probs = seq(0,0.01,0.001))),quantile(TDandFWH_Klet$FayandWusH_Klet,probs = seq(0,0.01,0.001)))
FWHPercKowa<-rbind(names(quantile(TDandFWH_Kowa$FayandWusH_Kowa,probs = seq(0,0.01,0.001))),quantile(TDandFWH_Kowa$FayandWusH_Kowa,probs = seq(0,0.01,0.001)))
TDPercKletall<-rbind(names(quantile(wholedata$TajimasD_Klet,probs = seq(0,0.01,0.001))),quantile(wholedata$TajimasD_Klet,probs = seq(0,0.01,0.001)))
TDPercKowaall<-rbind(names(quantile(wholedata$TajimasD_Kowa,probs = seq(0,0.01,0.001))),quantile(wholedata$TajimasD_Kowa,probs = seq(0,0.01,0.001)))
FWHPercKletall<-rbind(names(quantile(wholedata$FayandWusH_Klet,probs = seq(0,0.01,0.001))),quantile(wholedata$FayandWusH_Klet,probs = seq(0,0.01,0.001)))
FWHPercKowaall<-rbind(names(quantile(wholedata$FayandWusH_Kowa,probs = seq(0,0.01,0.001))),quantile(wholedata$FayandWusH_Kowa,probs = seq(0,0.01,0.001)))
DDPerc<-rbind(names(quantile(wholedata$DD,probs = seq(0,0.01,0.001))),quantile(wholedata$DD,probs = seq(0,0.01,0.001)))

Quantiles<-rbind(FSTPerc,AFDabsPerc,NIELSENPerc,DXYPerc,VARLDPerc,FLKPerc,TDPercKlet,TDPercKowa,FWHPercKlet,
FWHPercKowa,TDPercKletall,TDPercKowaall,FWHPercKletall,FWHPercKowaall,DDPerc)
Quantiles2<-cbind(c("FST","FST","AFDabs","AFDabs","NIELSEN","NIELSEN","DXY","DXY","VARLD","VARLD","FLK","FLK","TDKlet","TDKlet","TDKowa","TDKowa","FWHKlet"
,"FWHKlet","FWHKowa","FWHKowa","TDKletall","TDKletall","TDKowaall","TDKowaall","FWHKletall","FWHKletall","FWHKowaall","FWHKowaall","DD","DD"),Quantiles)

write.xlsx2(Quantiles2,"Quantiles_KletKowa_arenosa.xlsx",sheetName="Quantiles",col.names=FALSE,row.names=FALSE)



jpeg("Fst_Flk_KletKowa.jpeg", width=18, height=18, units="cm", res=1000)
plot(wholedata$Fst~wholedata$Flk)
abline(lm(wholedata$Fst~wholedata$Flk),col="red")
dev.off()

#save files in GeneslyrataKletKowa_1%outliers_TDmean_arenosa.xlsx as single txt file in Excel and load files
dataKletFST<-read.csv("Genes1percentFst_lyratagenome_TDmean_arenosa.txt",header=T,sep="\t")
dataKletNIELSEN<-read.csv("Genes1percentNielsen_lyratagenome_TDmean_arenosa.txt",header=T,sep="\t")
dataKletDXY<-read.csv("Genes1percentDxy_lyratagenome_TDmean_arenosa.txt",header=T,sep="\t")
dataKletFLK<-read.csv("Genes1percentFlk_lyratagenome_TDmean_arenosa.txt",header=T,sep="\t")
dataKletVARLD<-read.csv("Genes1percentVarLD_lyratagenome_TDmean_arenosa.txt",header=T,sep="\t")
dataKletAFDabs<-read.csv("Genes1percentAFD_lyratagenome_TDmean_arenosa.txt",header=T,sep="\t")
dataKletTD<-read.csv("Genes1percentTajimasD_Klet_lyratagenome_TDmean.txt",header=T,sep="\t")
dataKletDD<-read.csv("Genes1percentDD_lyratagenome_TDmean_arenosa.txt",header=T,sep="\t")

all<-unique(c(as.character(dataKletFST$Lyr_Gene),as.character(dataKletNIELSEN$Lyr_Gene),as.character(dataKletDXY$Lyr_Gene),as.character(dataKletFLK$Lyr_Gene),as.character(dataKletVARLD$Lyr_Gene),as.character(dataKletAFDabs$Lyr_Gene),as.character(dataKletTD$Lyr_Gene),as.character(dataKletDD$Lyr_Gene)))
write.table(all,"Lyrgenes_1perc_KK_arenosa.table",sep="\t",row.names=F)

#load SNPeff data
SNPeffdf<-read.table("SNPeff_1perc_KK_arenosa.table",header=T,sep="\t")

require(GenomicRanges)
SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]
KletFST_GRange<-GRanges(seqnames=dataKletFST$Chr,ranges=IRanges(start=dataKletFST$Window_start,end=dataKletFST$Window_end))
values(KletFST_GRange)<-dataKletFST[,1:30]
SNPeffdfMODIFIER_GRange<-GRanges(seqnames=SNPeffdfMODIFIER$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER$POS,end=SNPeffdfMODIFIER$POS))
values(SNPeffdfMODIFIER_GRange)<-SNPeffdfMODIFIER[,1:6]
FstallMODIFIER=mergeByOverlaps(KletFST_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
FstallMODIFIER_df=as.data.frame(FstallMODIFIER)
FstallMODIFIER_list<-FstallMODIFIER_df[as.character(FstallMODIFIER_df$Lyr_Gene)==as.character(FstallMODIFIER_df$Gene.1),]
dataKletFST$MODIFIER<-rep("NO",nrow(dataKletFST))
dataKletFST$MODIFIER<-ifelse(dataKletFST$Lyr_Gene%in%FstallMODIFIER_list$Lyr_Gene&dataKletFST$Lyr_Gene%in%FstallMODIFIER_list$Lyr_Gene,"MODIFIER",dataKletFST$MODIFIER)

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
KletFST_GRange<-GRanges(seqnames=dataKletFST$Chr,ranges=IRanges(start=dataKletFST$Window_start,end=dataKletFST$Window_end))
values(KletFST_GRange)<-dataKletFST[,1:30]
SNPeffdfMODERATE_GRange<-GRanges(seqnames=SNPeffdfMODERATE$CHROM,ranges=IRanges(start=SNPeffdfMODERATE$POS,end=SNPeffdfMODERATE$POS))
values(SNPeffdfMODERATE_GRange)<-SNPeffdfMODERATE[,1:6]
FstallMODERATE=mergeByOverlaps(KletFST_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
FstallMODERATE_df=as.data.frame(FstallMODERATE)
FstallMODERATE_list<-FstallMODERATE_df[as.character(FstallMODERATE_df$Lyr_Gene)==as.character(FstallMODERATE_df$Gene.1),]
dataKletFST$MODERATE<-rep("NO",nrow(dataKletFST))
dataKletFST$MODERATE<-ifelse(dataKletFST$Lyr_Gene%in%FstallMODERATE_list$Lyr_Gene&dataKletFST$Lyr_Gene%in%FstallMODERATE_list$Lyr_Gene,"MODERATE",dataKletFST$MODERATE)

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
KletFST_GRange<-GRanges(seqnames=dataKletFST$Chr,ranges=IRanges(start=dataKletFST$Window_start,end=dataKletFST$Window_end))
values(KletFST_GRange)<-dataKletFST[,1:30]
SNPeffdfHIGH_GRange<-GRanges(seqnames=SNPeffdfHIGH$CHROM,ranges=IRanges(start=SNPeffdfHIGH$POS,end=SNPeffdfHIGH$POS))
values(SNPeffdfHIGH_GRange)<-SNPeffdfHIGH[,1:6]
FstallHIGH=mergeByOverlaps(KletFST_GRange,SNPeffdfHIGH_GRange,type=c("any"))
FstallHIGH_df=as.data.frame(FstallHIGH)
FstallHIGH_list<-FstallHIGH_df[as.character(FstallHIGH_df$Lyr_Gene)==as.character(FstallHIGH_df$Gene.1),]
dataKletFST$HIGH<-rep("NO",nrow(dataKletFST))
dataKletFST$HIGH<-ifelse(dataKletFST$Lyr_Gene%in%FstallHIGH_list$Lyr_Gene&dataKletFST$Lyr_Gene%in%FstallHIGH_list$Lyr_Gene,"HIGH",dataKletFST$HIGH)

write.table(dataKletFST,"Summary_genes_1percent_KletKowa_FST_arenosa.csv",col.names=TRUE,row.names=FALSE,sep=";")

SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]
KletNIELSEN_GRange<-GRanges(seqnames=dataKletNIELSEN$Chr,ranges=IRanges(start=dataKletNIELSEN$Window_start,end=dataKletNIELSEN$Window_end))
values(KletNIELSEN_GRange)<-dataKletNIELSEN[,1:30]
SNPeffdfMODIFIER_GRange<-GRanges(seqnames=SNPeffdfMODIFIER$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER$POS,end=SNPeffdfMODIFIER$POS))
values(SNPeffdfMODIFIER_GRange)<-SNPeffdfMODIFIER[,1:6]
NielsenallMODIFIER=mergeByOverlaps(KletNIELSEN_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
NielsenallMODIFIER_df=as.data.frame(NielsenallMODIFIER)
NielsenallMODIFIER_list<-NielsenallMODIFIER_df[as.character(NielsenallMODIFIER_df$Lyr_Gene)==as.character(NielsenallMODIFIER_df$Gene.1),]
dataKletNIELSEN$MODIFIER<-rep("NO",nrow(dataKletNIELSEN))
dataKletNIELSEN$MODIFIER<-ifelse(dataKletNIELSEN$Lyr_Gene%in%NielsenallMODIFIER_list$Lyr_Gene&dataKletNIELSEN$Lyr_Gene%in%NielsenallMODIFIER_list$Lyr_Gene,"MODIFIER",dataKletNIELSEN$MODIFIER)

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
KletNIELSEN_GRange<-GRanges(seqnames=dataKletNIELSEN$Chr,ranges=IRanges(start=dataKletNIELSEN$Window_start,end=dataKletNIELSEN$Window_end))
values(KletNIELSEN_GRange)<-dataKletNIELSEN[,1:30]
SNPeffdfMODERATE_GRange<-GRanges(seqnames=SNPeffdfMODERATE$CHROM,ranges=IRanges(start=SNPeffdfMODERATE$POS,end=SNPeffdfMODERATE$POS))
values(SNPeffdfMODERATE_GRange)<-SNPeffdfMODERATE[,1:6]
NielsenallMODERATE=mergeByOverlaps(KletNIELSEN_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
NielsenallMODERATE_df=as.data.frame(NielsenallMODERATE)
NielsenallMODERATE_list<-NielsenallMODERATE_df[as.character(NielsenallMODERATE_df$Lyr_Gene)==as.character(NielsenallMODERATE_df$Gene.1),]
dataKletNIELSEN$MODERATE<-rep("NO",nrow(dataKletNIELSEN))
dataKletNIELSEN$MODERATE<-ifelse(dataKletNIELSEN$Lyr_Gene%in%NielsenallMODERATE_list$Lyr_Gene&dataKletNIELSEN$Lyr_Gene%in%NielsenallMODERATE_list$Lyr_Gene,"MODERATE",dataKletNIELSEN$MODERATE)

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
KletNIELSEN_GRange<-GRanges(seqnames=dataKletNIELSEN$Chr,ranges=IRanges(start=dataKletNIELSEN$Window_start,end=dataKletNIELSEN$Window_end))
values(KletNIELSEN_GRange)<-dataKletNIELSEN[,1:30]
SNPeffdfHIGH_GRange<-GRanges(seqnames=SNPeffdfHIGH$CHROM,ranges=IRanges(start=SNPeffdfHIGH$POS,end=SNPeffdfHIGH$POS))
values(SNPeffdfHIGH_GRange)<-SNPeffdfHIGH[,1:6]
NielsenallHIGH=mergeByOverlaps(KletNIELSEN_GRange,SNPeffdfHIGH_GRange,type=c("any"))
NielsenallHIGH_df=as.data.frame(NielsenallHIGH)
NielsenallHIGH_list<-NielsenallHIGH_df[as.character(NielsenallHIGH_df$Lyr_Gene)==as.character(NielsenallHIGH_df$Gene.1),]
dataKletNIELSEN$HIGH<-rep("NO",nrow(dataKletNIELSEN))
dataKletNIELSEN$HIGH<-ifelse(dataKletNIELSEN$Lyr_Gene%in%NielsenallHIGH_list$Lyr_Gene&dataKletNIELSEN$Lyr_Gene%in%NielsenallHIGH_list$Lyr_Gene,"HIGH",dataKletNIELSEN$HIGH)

write.table(dataKletNIELSEN,"Summary_genes_1percent_KletKowa_NIELSEN_arenosa.csv",col.names=TRUE,row.names=FALSE,sep=";")


SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]
KletDD_GRange<-GRanges(seqnames=dataKletDD$Chr,ranges=IRanges(start=dataKletDD$Window_start,end=dataKletDD$Window_end))
values(KletDD_GRange)<-dataKletDD[,1:30]
SNPeffdfMODIFIER_GRange<-GRanges(seqnames=SNPeffdfMODIFIER$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER$POS,end=SNPeffdfMODIFIER$POS))
values(SNPeffdfMODIFIER_GRange)<-SNPeffdfMODIFIER[,1:6]
DDallMODIFIER=mergeByOverlaps(KletDD_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
DDallMODIFIER_df=as.data.frame(DDallMODIFIER)
DDallMODIFIER_list<-DDallMODIFIER_df[as.character(DDallMODIFIER_df$Lyr_Gene)==as.character(DDallMODIFIER_df$Gene.1),]
dataKletDD$MODIFIER<-rep("NO",nrow(dataKletDD))
dataKletDD$MODIFIER<-ifelse(dataKletDD$Lyr_Gene%in%DDallMODIFIER_list$Lyr_Gene&dataKletDD$Lyr_Gene%in%DDallMODIFIER_list$Lyr_Gene,"MODIFIER",dataKletDD$MODIFIER)

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
KletDD_GRange<-GRanges(seqnames=dataKletDD$Chr,ranges=IRanges(start=dataKletDD$Window_start,end=dataKletDD$Window_end))
values(KletDD_GRange)<-dataKletDD[,1:30]
SNPeffdfMODERATE_GRange<-GRanges(seqnames=SNPeffdfMODERATE$CHROM,ranges=IRanges(start=SNPeffdfMODERATE$POS,end=SNPeffdfMODERATE$POS))
values(SNPeffdfMODERATE_GRange)<-SNPeffdfMODERATE[,1:6]
DDallMODERATE=mergeByOverlaps(KletDD_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
DDallMODERATE_df=as.data.frame(DDallMODERATE)
DDallMODERATE_list<-DDallMODERATE_df[as.character(DDallMODERATE_df$Lyr_Gene)==as.character(DDallMODERATE_df$Gene.1),]
dataKletDD$MODERATE<-rep("NO",nrow(dataKletDD))
dataKletDD$MODERATE<-ifelse(dataKletDD$Lyr_Gene%in%DDallMODERATE_list$Lyr_Gene&dataKletDD$Lyr_Gene%in%DDallMODERATE_list$Lyr_Gene,"MODERATE",dataKletDD$MODERATE)

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
KletDD_GRange<-GRanges(seqnames=dataKletDD$Chr,ranges=IRanges(start=dataKletDD$Window_start,end=dataKletDD$Window_end))
values(KletDD_GRange)<-dataKletDD[,1:30]
SNPeffdfHIGH_GRange<-GRanges(seqnames=SNPeffdfHIGH$CHROM,ranges=IRanges(start=SNPeffdfHIGH$POS,end=SNPeffdfHIGH$POS))
values(SNPeffdfHIGH_GRange)<-SNPeffdfHIGH[,1:6]
DDallHIGH=mergeByOverlaps(KletDD_GRange,SNPeffdfHIGH_GRange,type=c("any"))
DDallHIGH_df=as.data.frame(DDallHIGH)
DDallHIGH_list<-DDallHIGH_df[as.character(DDallHIGH_df$Lyr_Gene)==as.character(DDallHIGH_df$Gene.1),]
dataKletDD$HIGH<-rep("NO",nrow(dataKletDD))
dataKletDD$HIGH<-ifelse(dataKletDD$Lyr_Gene%in%DDallHIGH_list$Lyr_Gene&dataKletDD$Lyr_Gene%in%DDallHIGH_list$Lyr_Gene,"HIGH",dataKletDD$HIGH)

write.table(dataKletDD,"Summary_genes_1percent_KletKowa_DD_arenosa.csv",col.names=TRUE,row.names=FALSE,sep=";")

SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]
KletDXY_GRange<-GRanges(seqnames=dataKletDXY$Chr,ranges=IRanges(start=dataKletDXY$Window_start,end=dataKletDXY$Window_end))
values(KletDXY_GRange)<-dataKletDXY[,1:30]
SNPeffdfMODIFIER_GRange<-GRanges(seqnames=SNPeffdfMODIFIER$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER$POS,end=SNPeffdfMODIFIER$POS))
values(SNPeffdfMODIFIER_GRange)<-SNPeffdfMODIFIER[,1:6]
DxyallMODIFIER=mergeByOverlaps(KletDXY_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
DxyallMODIFIER_df=as.data.frame(DxyallMODIFIER)
DxyallMODIFIER_list<-DxyallMODIFIER_df[as.character(DxyallMODIFIER_df$Lyr_Gene)==as.character(DxyallMODIFIER_df$Gene.1),]
dataKletDXY$MODIFIER<-rep("NO",nrow(dataKletDXY))
dataKletDXY$MODIFIER<-ifelse(dataKletDXY$Lyr_Gene%in%DxyallMODIFIER_list$Lyr_Gene&dataKletDXY$Lyr_Gene%in%DxyallMODIFIER_list$Lyr_Gene,"MODIFIER",dataKletDXY$MODIFIER)

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
KletDXY_GRange<-GRanges(seqnames=dataKletDXY$Chr,ranges=IRanges(start=dataKletDXY$Window_start,end=dataKletDXY$Window_end))
values(KletDXY_GRange)<-dataKletDXY[,1:30]
SNPeffdfMODERATE_GRange<-GRanges(seqnames=SNPeffdfMODERATE$CHROM,ranges=IRanges(start=SNPeffdfMODERATE$POS,end=SNPeffdfMODERATE$POS))
values(SNPeffdfMODERATE_GRange)<-SNPeffdfMODERATE[,1:6]
DxyallMODERATE=mergeByOverlaps(KletDXY_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
DxyallMODERATE_df=as.data.frame(DxyallMODERATE)
DxyallMODERATE_list<-DxyallMODERATE_df[as.character(DxyallMODERATE_df$Lyr_Gene)==as.character(DxyallMODERATE_df$Gene.1),]
dataKletDXY$MODERATE<-rep("NO",nrow(dataKletDXY))
dataKletDXY$MODERATE<-ifelse(dataKletDXY$Lyr_Gene%in%DxyallMODERATE_list$Lyr_Gene&dataKletDXY$Lyr_Gene%in%DxyallMODERATE_list$Lyr_Gene,"MODERATE",dataKletDXY$MODERATE)

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
KletDXY_GRange<-GRanges(seqnames=dataKletDXY$Chr,ranges=IRanges(start=dataKletDXY$Window_start,end=dataKletDXY$Window_end))
values(KletDXY_GRange)<-dataKletDXY[,1:30]
SNPeffdfHIGH_GRange<-GRanges(seqnames=SNPeffdfHIGH$CHROM,ranges=IRanges(start=SNPeffdfHIGH$POS,end=SNPeffdfHIGH$POS))
values(SNPeffdfHIGH_GRange)<-SNPeffdfHIGH[,1:6]
DxyallHIGH=mergeByOverlaps(KletDXY_GRange,SNPeffdfHIGH_GRange,type=c("any"))
DxyallHIGH_df=as.data.frame(DxyallHIGH)
DxyallHIGH_list<-DxyallHIGH_df[as.character(DxyallHIGH_df$Lyr_Gene)==as.character(DxyallHIGH_df$Gene.1),]
dataKletDXY$HIGH<-rep("NO",nrow(dataKletDXY))
dataKletDXY$HIGH<-ifelse(dataKletDXY$Lyr_Gene%in%DxyallHIGH_list$Lyr_Gene&dataKletDXY$Lyr_Gene%in%DxyallHIGH_list$Lyr_Gene,"HIGH",dataKletDXY$HIGH)

write.table(dataKletDXY,"Summary_genes_1percent_KletKowa_DXY_arenosa.csv",col.names=TRUE,row.names=FALSE,sep=";")

SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]
KletFLK_GRange<-GRanges(seqnames=dataKletFLK$Chr,ranges=IRanges(start=dataKletFLK$Window_start,end=dataKletFLK$Window_end))
values(KletFLK_GRange)<-dataKletFLK[,1:30]
SNPeffdfMODIFIER_GRange<-GRanges(seqnames=SNPeffdfMODIFIER$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER$POS,end=SNPeffdfMODIFIER$POS))
values(SNPeffdfMODIFIER_GRange)<-SNPeffdfMODIFIER[,1:6]
FlkallMODIFIER=mergeByOverlaps(KletFLK_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
FlkallMODIFIER_df=as.data.frame(FlkallMODIFIER)
FlkallMODIFIER_list<-FlkallMODIFIER_df[as.character(FlkallMODIFIER_df$Lyr_Gene)==as.character(FlkallMODIFIER_df$Gene.1),]
dataKletFLK$MODIFIER<-rep("NO",nrow(dataKletFLK))
dataKletFLK$MODIFIER<-ifelse(dataKletFLK$Lyr_Gene%in%FlkallMODIFIER_list$Lyr_Gene&dataKletFLK$Lyr_Gene%in%FlkallMODIFIER_list$Lyr_Gene,"MODIFIER",dataKletFLK$MODIFIER)

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
KletFLK_GRange<-GRanges(seqnames=dataKletFLK$Chr,ranges=IRanges(start=dataKletFLK$Window_start,end=dataKletFLK$Window_end))
values(KletFLK_GRange)<-dataKletFLK[,1:30]
SNPeffdfMODERATE_GRange<-GRanges(seqnames=SNPeffdfMODERATE$CHROM,ranges=IRanges(start=SNPeffdfMODERATE$POS,end=SNPeffdfMODERATE$POS))
values(SNPeffdfMODERATE_GRange)<-SNPeffdfMODERATE[,1:6]
FlkallMODERATE=mergeByOverlaps(KletFLK_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
FlkallMODERATE_df=as.data.frame(FlkallMODERATE)
FlkallMODERATE_list<-FlkallMODERATE_df[as.character(FlkallMODERATE_df$Lyr_Gene)==as.character(FlkallMODERATE_df$Gene.1),]
dataKletFLK$MODERATE<-rep("NO",nrow(dataKletFLK))
dataKletFLK$MODERATE<-ifelse(dataKletFLK$Lyr_Gene%in%FlkallMODERATE_list$Lyr_Gene&dataKletFLK$Lyr_Gene%in%FlkallMODERATE_list$Lyr_Gene,"MODERATE",dataKletFLK$MODERATE)

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
KletFLK_GRange<-GRanges(seqnames=dataKletFLK$Chr,ranges=IRanges(start=dataKletFLK$Window_start,end=dataKletFLK$Window_end))
values(KletFLK_GRange)<-dataKletFLK[,1:30]
SNPeffdfHIGH_GRange<-GRanges(seqnames=SNPeffdfHIGH$CHROM,ranges=IRanges(start=SNPeffdfHIGH$POS,end=SNPeffdfHIGH$POS))
values(SNPeffdfHIGH_GRange)<-SNPeffdfHIGH[,1:6]
FlkallHIGH=mergeByOverlaps(KletFLK_GRange,SNPeffdfHIGH_GRange,type=c("any"))
FlkallHIGH_df=as.data.frame(FlkallHIGH)
FlkallHIGH_list<-FlkallHIGH_df[as.character(FlkallHIGH_df$Lyr_Gene)==as.character(FlkallHIGH_df$Gene.1),]
dataKletFLK$HIGH<-rep("NO",nrow(dataKletFLK))
dataKletFLK$HIGH<-ifelse(dataKletFLK$Lyr_Gene%in%FlkallHIGH_list$Lyr_Gene&dataKletFLK$Lyr_Gene%in%FlkallHIGH_list$Lyr_Gene,"HIGH",dataKletFLK$HIGH)

write.table(dataKletFLK,"Summary_genes_1percent_KletKowa_FLK_arenosa.csv",col.names=TRUE,row.names=FALSE,sep=";")

SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]
KletAFDabs_GRange<-GRanges(seqnames=dataKletAFDabs$Chr,ranges=IRanges(start=dataKletAFDabs$Window_start,end=dataKletAFDabs$Window_end))
values(KletAFDabs_GRange)<-dataKletAFDabs[,1:30]
SNPeffdfMODIFIER_GRange<-GRanges(seqnames=SNPeffdfMODIFIER$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER$POS,end=SNPeffdfMODIFIER$POS))
values(SNPeffdfMODIFIER_GRange)<-SNPeffdfMODIFIER[,1:6]
AFDabsallMODIFIER=mergeByOverlaps(KletAFDabs_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
AFDabsallMODIFIER_df=as.data.frame(AFDabsallMODIFIER)
AFDabsallMODIFIER_list<-AFDabsallMODIFIER_df[as.character(AFDabsallMODIFIER_df$Lyr_Gene)==as.character(AFDabsallMODIFIER_df$Gene.1),]
dataKletAFDabs$MODIFIER<-rep("NO",nrow(dataKletAFDabs))
dataKletAFDabs$MODIFIER<-ifelse(dataKletAFDabs$Lyr_Gene%in%AFDabsallMODIFIER_list$Lyr_Gene&dataKletAFDabs$Lyr_Gene%in%AFDabsallMODIFIER_list$Lyr_Gene,"MODIFIER",dataKletAFDabs$MODIFIER)

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
KletAFDabs_GRange<-GRanges(seqnames=dataKletAFDabs$Chr,ranges=IRanges(start=dataKletAFDabs$Window_start,end=dataKletAFDabs$Window_end))
values(KletAFDabs_GRange)<-dataKletAFDabs[,1:30]
SNPeffdfMODERATE_GRange<-GRanges(seqnames=SNPeffdfMODERATE$CHROM,ranges=IRanges(start=SNPeffdfMODERATE$POS,end=SNPeffdfMODERATE$POS))
values(SNPeffdfMODERATE_GRange)<-SNPeffdfMODERATE[,1:6]
AFDabsallMODERATE=mergeByOverlaps(KletAFDabs_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
AFDabsallMODERATE_df=as.data.frame(AFDabsallMODERATE)
AFDabsallMODERATE_list<-AFDabsallMODERATE_df[as.character(AFDabsallMODERATE_df$Lyr_Gene)==as.character(AFDabsallMODERATE_df$Gene.1),]
dataKletAFDabs$MODERATE<-rep("NO",nrow(dataKletAFDabs))
dataKletAFDabs$MODERATE<-ifelse(dataKletAFDabs$Lyr_Gene%in%AFDabsallMODERATE_list$Lyr_Gene&dataKletAFDabs$Lyr_Gene%in%AFDabsallMODERATE_list$Lyr_Gene,"MODERATE",dataKletAFDabs$MODERATE)

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
KletAFDabs_GRange<-GRanges(seqnames=dataKletAFDabs$Chr,ranges=IRanges(start=dataKletAFDabs$Window_start,end=dataKletAFDabs$Window_end))
values(KletAFDabs_GRange)<-dataKletAFDabs[,1:30]
SNPeffdfHIGH_GRange<-GRanges(seqnames=SNPeffdfHIGH$CHROM,ranges=IRanges(start=SNPeffdfHIGH$POS,end=SNPeffdfHIGH$POS))
values(SNPeffdfHIGH_GRange)<-SNPeffdfHIGH[,1:6]
AFDabsallHIGH=mergeByOverlaps(KletAFDabs_GRange,SNPeffdfHIGH_GRange,type=c("any"))
AFDabsallHIGH_df=as.data.frame(AFDabsallHIGH)
AFDabsallHIGH_list<-AFDabsallHIGH_df[as.character(AFDabsallHIGH_df$Lyr_Gene)==as.character(AFDabsallHIGH_df$Gene.1),]
dataKletAFDabs$HIGH<-rep("NO",nrow(dataKletAFDabs))
dataKletAFDabs$HIGH<-ifelse(dataKletAFDabs$Lyr_Gene%in%AFDabsallHIGH_list$Lyr_Gene&dataKletAFDabs$Lyr_Gene%in%AFDabsallHIGH_list$Lyr_Gene,"HIGH",dataKletAFDabs$HIGH)

write.table(dataKletAFDabs,"Summary_genes_1percent_KletKowa_AFDabs_arenosa.csv",col.names=TRUE,row.names=FALSE,sep=";")


SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]
KletVARLD_GRange<-GRanges(seqnames=dataKletVARLD$Chr,ranges=IRanges(start=dataKletVARLD$Window_start,end=dataKletVARLD$Window_end))
values(KletVARLD_GRange)<-dataKletVARLD[,1:30]
SNPeffdfMODIFIER_GRange<-GRanges(seqnames=SNPeffdfMODIFIER$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER$POS,end=SNPeffdfMODIFIER$POS))
values(SNPeffdfMODIFIER_GRange)<-SNPeffdfMODIFIER[,1:6]
VARLDallMODIFIER=mergeByOverlaps(KletVARLD_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
VARLDallMODIFIER_df=as.data.frame(VARLDallMODIFIER)
VARLDallMODIFIER_list<-VARLDallMODIFIER_df[as.character(VARLDallMODIFIER_df$Lyr_Gene)==as.character(VARLDallMODIFIER_df$Gene.1),]
dataKletVARLD$MODIFIER<-rep("NO",nrow(dataKletVARLD))
dataKletVARLD$MODIFIER<-ifelse(dataKletVARLD$Lyr_Gene%in%VARLDallMODIFIER_list$Lyr_Gene&dataKletVARLD$Lyr_Gene%in%VARLDallMODIFIER_list$Lyr_Gene,"MODIFIER",dataKletVARLD$MODIFIER)

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
KletVARLD_GRange<-GRanges(seqnames=dataKletVARLD$Chr,ranges=IRanges(start=dataKletVARLD$Window_start,end=dataKletVARLD$Window_end))
values(KletVARLD_GRange)<-dataKletVARLD[,1:30]
SNPeffdfMODERATE_GRange<-GRanges(seqnames=SNPeffdfMODERATE$CHROM,ranges=IRanges(start=SNPeffdfMODERATE$POS,end=SNPeffdfMODERATE$POS))
values(SNPeffdfMODERATE_GRange)<-SNPeffdfMODERATE[,1:6]
VARLDallMODERATE=mergeByOverlaps(KletVARLD_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
VARLDallMODERATE_df=as.data.frame(VARLDallMODERATE)
VARLDallMODERATE_list<-VARLDallMODERATE_df[as.character(VARLDallMODERATE_df$Lyr_Gene)==as.character(VARLDallMODERATE_df$Gene.1),]
dataKletVARLD$MODERATE<-rep("NO",nrow(dataKletVARLD))
dataKletVARLD$MODERATE<-ifelse(dataKletVARLD$Lyr_Gene%in%VARLDallMODERATE_list$Lyr_Gene&dataKletVARLD$Lyr_Gene%in%VARLDallMODERATE_list$Lyr_Gene,"MODERATE",dataKletVARLD$MODERATE)

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
KletVARLD_GRange<-GRanges(seqnames=dataKletVARLD$Chr,ranges=IRanges(start=dataKletVARLD$Window_start,end=dataKletVARLD$Window_end))
values(KletVARLD_GRange)<-dataKletVARLD[,1:30]
SNPeffdfHIGH_GRange<-GRanges(seqnames=SNPeffdfHIGH$CHROM,ranges=IRanges(start=SNPeffdfHIGH$POS,end=SNPeffdfHIGH$POS))
values(SNPeffdfHIGH_GRange)<-SNPeffdfHIGH[,1:6]
VARLDallHIGH=mergeByOverlaps(KletVARLD_GRange,SNPeffdfHIGH_GRange,type=c("any"))
VARLDallHIGH_df=as.data.frame(VARLDallHIGH)
VARLDallHIGH_list<-VARLDallHIGH_df[as.character(VARLDallHIGH_df$Lyr_Gene)==as.character(VARLDallHIGH_df$Gene.1),]
dataKletVARLD$HIGH<-rep("NO",nrow(dataKletVARLD))
dataKletVARLD$HIGH<-ifelse(dataKletVARLD$Lyr_Gene%in%VARLDallHIGH_list$Lyr_Gene&dataKletVARLD$Lyr_Gene%in%VARLDallHIGH_list$Lyr_Gene,"HIGH",dataKletVARLD$HIGH)

write.table(dataKletVARLD,"Summary_genes_1percent_KletKowa_VARLD_arenosa.csv",col.names=TRUE,row.names=FALSE,sep=";")

SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]
KletTD_GRange<-GRanges(seqnames=dataKletTD$Chr,ranges=IRanges(start=dataKletTD$Window_start,end=dataKletTD$Window_end))
values(KletTD_GRange)<-dataKletTD[,1:19]
SNPeffdfMODIFIER_GRange<-GRanges(seqnames=SNPeffdfMODIFIER$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER$POS,end=SNPeffdfMODIFIER$POS))
values(SNPeffdfMODIFIER_GRange)<-SNPeffdfMODIFIER[,1:6]
TDallMODIFIER=mergeByOverlaps(KletTD_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
TDallMODIFIER_df=as.data.frame(TDallMODIFIER)
TDallMODIFIER_list<-TDallMODIFIER_df[as.character(TDallMODIFIER_df$Lyr_Gene)==as.character(TDallMODIFIER_df$Gene.1),]
dataKletTD$MODIFIER<-rep("NO",nrow(dataKletTD))
dataKletTD$MODIFIER<-ifelse(dataKletTD$Lyr_Gene%in%TDallMODIFIER_list$Lyr_Gene&dataKletTD$Lyr_Gene%in%TDallMODIFIER_list$Lyr_Gene,"MODIFIER",dataKletTD$MODIFIER)

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
KletTD_GRange<-GRanges(seqnames=dataKletTD$Chr,ranges=IRanges(start=dataKletTD$Window_start,end=dataKletTD$Window_end))
values(KletTD_GRange)<-dataKletTD[,1:19]
SNPeffdfMODERATE_GRange<-GRanges(seqnames=SNPeffdfMODERATE$CHROM,ranges=IRanges(start=SNPeffdfMODERATE$POS,end=SNPeffdfMODERATE$POS))
values(SNPeffdfMODERATE_GRange)<-SNPeffdfMODERATE[,1:6]
TDallMODERATE=mergeByOverlaps(KletTD_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
TDallMODERATE_df=as.data.frame(TDallMODERATE)
TDallMODERATE_list<-TDallMODERATE_df[as.character(TDallMODERATE_df$Lyr_Gene)==as.character(TDallMODERATE_df$Gene.1),]
dataKletTD$MODERATE<-rep("NO",nrow(dataKletTD))
dataKletTD$MODERATE<-ifelse(dataKletTD$Lyr_Gene%in%TDallMODERATE_list$Lyr_Gene&dataKletTD$Lyr_Gene%in%TDallMODERATE_list$Lyr_Gene,"MODERATE",dataKletTD$MODERATE)

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
KletTD_GRange<-GRanges(seqnames=dataKletTD$Chr,ranges=IRanges(start=dataKletTD$Window_start,end=dataKletTD$Window_end))
values(KletTD_GRange)<-dataKletTD[,1:19]
SNPeffdfHIGH_GRange<-GRanges(seqnames=SNPeffdfHIGH$CHROM,ranges=IRanges(start=SNPeffdfHIGH$POS,end=SNPeffdfHIGH$POS))
values(SNPeffdfHIGH_GRange)<-SNPeffdfHIGH[,1:6]
TDallHIGH=mergeByOverlaps(KletTD_GRange,SNPeffdfHIGH_GRange,type=c("any"))
TDallHIGH_df=as.data.frame(TDallHIGH)
TDallHIGH_list<-TDallHIGH_df[as.character(TDallHIGH_df$Lyr_Gene)==as.character(TDallHIGH_df$Gene.1),]
dataKletTD$HIGH<-rep("NO",nrow(dataKletTD))
dataKletTD$HIGH<-ifelse(dataKletTD$Lyr_Gene%in%TDallHIGH_list$Lyr_Gene&dataKletTD$Lyr_Gene%in%TDallHIGH_list$Lyr_Gene,"HIGH",dataKletTD$HIGH)

write.table(dataKletTD,"Summary_genes_1percent_KletKowa_TD_arenosa.csv",col.names=TRUE,row.names=FALSE,sep=";")

require(xlsx)
#values taken from Quantiles_KletKowa_arenosa.xlsx

FST<-read.csv("Summary_genes_1percent_KletKowa_FST_arenosa.csv",header=T,sep=";")
FST01<-FST[FST$Fst>=0.159584258881707,]
write.xlsx2(FST01,"Genes_01percent_KletKowa_arenosa.xlsx",sheetName="Fst_0_1%",col.names=TRUE,row.names=FALSE)

DXY<-read.csv("Summary_genes_1percent_KletKowa_DXY_arenosa.csv",header=T,sep=";")
DXY01<-DXY[DXY$Dxy>=0.419311666666666,]
write.xlsx2(DXY01,"Genes_01percent_KletKowa_arenosa.xlsx",sheetName="Dxy_0_1%",col.names=TRUE,row.names=FALSE,append=TRUE)

AFDabs<-read.csv("Summary_genes_1percent_KletKowa_AFDabs_arenosa.csv",header=T,sep=";")
AFDabs01<-AFDabs[AFDabs$AFDabs>=0.288629444444442,]
write.xlsx2(AFDabs01,"Genes_01percent_KletKowa_arenosa.xlsx",sheetName="AFDabs_0_1%",col.names=TRUE,row.names=FALSE,append=TRUE)

NIELSEN<-read.csv("Summary_genes_1percent_KletKowa_NIELSEN_arenosa.csv",header=T,sep=";")
NIELSEN01<-NIELSEN[NIELSEN$Nielsen>=243.271265489257,]
write.xlsx2(NIELSEN01,"Genes_01percent_KletKowa_arenosa.xlsx",sheetName="Nielsen_0_1%",col.names=TRUE,row.names=FALSE,append=TRUE)

VARLD<-read.csv("Summary_genes_1percent_KletKowa_VARLD_arenosa.csv",header=T,sep=";")
VARLD01<-VARLD[VARLD$VarLD>=34.4780176,]
write.xlsx2(VARLD01,"Genes_01percent_KletKowa_arenosa.xlsx",sheetName="VarLD_0_1%",col.names=TRUE,row.names=FALSE,append=TRUE)

FLK<-read.csv("Summary_genes_1percent_KletKowa_FLK_arenosa.csv",header=T,sep=";")
FLK01<-FLK[FLK$Flk>=3.29076116503269&FLK$Flk_pvalue<0.05,]
write.xlsx2(FLK01,"Genes_01percent_KletKowa_arenosa.xlsx",sheetName="Flk_0_1%",col.names=TRUE,row.names=FALSE,append=TRUE)

DD<-read.csv("Summary_genes_1percent_KletKowa_DD_arenosa.csv",header=T,sep=";")
DD01<-DD[DD$DD<=-0.147369851025826,]
write.xlsx2(DD01,"Genes_01percent_KletKowa_arenosa.xlsx",sheetName="DD_0_1%",col.names=TRUE,row.names=FALSE,append=TRUE)

#create and annotate allele frequency and metrics plots 
#check if different plot rating between metrics

FST<-read.xlsx2("Genes_01percent_KletKowa.xlsx",1)
DXY<-read.xlsx2("Genes_01percent_KletKowa.xlsx",2)
AFDabs<-read.xlsx2("Genes_01percent_KletKowa.xlsx",3)
NIELSEN<-read.xlsx2("Genes_01percent_KletKowa.xlsx",4)
VARLD<-read.xlsx2("Genes_01percent_KletKowa.xlsx",5)
FLK<-read.xlsx2("Genes_01percent_KletKowa.xlsx",6)
DD<-read.xlsx2("Genes_01percent_KletKowa.xlsx",7)

FSTDXY<-merge(FST,DXY,by="Gene")
nomatchFSTDXY<-ifelse(FSTDXY$AF_Plot.x!=FSTDXY$AF_Plot.y,as.character(FSTDXY$Gene),NA)

AFDabsDXY<-merge(AFDabs,DXY,by="Gene")
nomatchAFDabsDXY<-ifelse(AFDabsDXY$AF_Plot.x!=AFDabsDXY$AF_Plot.y,as.character(AFDabsDXY$Gene),NA)

AFDabsNIELSEN<-merge(AFDabs,NIELSEN,by="Gene")
nomatchAFDabsNIELSEN<-ifelse(AFDabsNIELSEN$AF_Plot.x!=AFDabsNIELSEN$AF_Plot.y,as.character(AFDabsNIELSEN$Gene),NA)

AFDabsVARLD<-merge(AFDabs,VARLD,by="Gene")
nomatchAFDabsVARLD<-ifelse(as.character(AFDabsVARLD$AF_Plot.x)!=as.character(AFDabsVARLD$AF_Plot.y),as.character(AFDabsVARLD$Gene),NA)

AFDabsFLK<-merge(AFDabs,FLK,by="Gene")
nomatchAFDabsFLK<-ifelse(AFDabsFLK$AF_Plot.x!=AFDabsFLK$AF_Plot.y,as.character(AFDabsFLK$Gene),NA)

AFDabsDD<-merge(AFDabs,DD,by="Gene")
nomatchAFDabsDD<-ifelse(as.character(AFDabsDD$AF_Plot.x)!=as.character(AFDabsDD$AF_Plot.y),as.character(AFDabsDD$Gene),NA)

all<-Reduce(function(x,y) merge(x,y,all=TRUE,by="Gene"),list(FST,DXY,AFDabs,NIELSEN,VARLD,FLK))
nomatchall<-ifelse(all$AF_Plot.x!=all$AF_Plot.y,as.character(all$Gene),NA)


FST<-read.xlsx2("Genes_01percent_KletKowa.xlsx",1)
DXY<-read.xlsx2("Genes_01percent_KletKowa.xlsx",2)
AFDabs<-read.xlsx2("Genes_01percent_KletKowa.xlsx",3)
NIELSEN<-read.xlsx2("Genes_01percent_KletKowa.xlsx",4)
VARLD<-read.xlsx2("Genes_01percent_KletKowa.xlsx",5)
FLK<-read.xlsx2("Genes_01percent_KletKowa.xlsx",6)
DD<-read.xlsx2("Genes_01percent_KletKowa.xlsx",7)

FSTDXY<-merge(FST,DXY,by="Gene")
nomatchFSTDXY<-ifelse(FSTDXY$Metrics_Plot.x!=FSTDXY$Metrics_Plot.y,as.character(FSTDXY$Gene),NA)

AFDabsDXY<-merge(AFDabs,DXY,by="Gene")
nomatchAFDabsDXY<-ifelse(AFDabsDXY$Metrics_Plot.x!=AFDabsDXY$Metrics_Plot.y,as.character(AFDabsDXY$Gene),NA)

AFDabsNIELSEN<-merge(AFDabs,NIELSEN,by="Gene")
nomatchAFDabsNIELSEN<-ifelse(AFDabsNIELSEN$Metrics_Plot.x!=AFDabsNIELSEN$Metrics_Plot.y,as.character(AFDabsNIELSEN$Gene),NA)

AFDabsVARLD<-merge(AFDabs,VARLD,by="Gene")
nomatchAFDabsVARLD<-ifelse(as.character(AFDabsVARLD$Metrics_Plot.x)!=as.character(AFDabsVARLD$Metrics_Plot.y),as.character(AFDabsVARLD$Gene),NA)

AFDabsFLK<-merge(AFDabs,FLK,by="Gene")
nomatchAFDabsFLK<-ifelse(AFDabsFLK$Metrics_Plot.x!=AFDabsFLK$Metrics_Plot.y,as.character(AFDabsFLK$Gene),NA)

AFDabsDD<-merge(AFDabs,DD,by="Gene")
nomatchAFDabsDD<-ifelse(as.character(AFDabsDD$Metrics_Plot.x)!=as.character(AFDabsDD$Metrics_Plot.y),as.character(AFDabsDD$Gene),NA)

all<-Reduce(function(x,y) merge(x,y,all=TRUE,by="Gene"),list(FST,DXY,AFDabs,NIELSEN,VARLD,FLK))
nomatchall<-ifelse(all$Metrics_Plot.x!=all$Metrics_Plot.y,as.character(all$Gene),NA)


FSTred<-FST[!(FST$AF_Plot=="D"|FST$Metrics_Plot=="C"),] 
DXYred<-DXY[!(DXY$AF_Plot=="D"|DXY$Metrics_Plot=="C"),]
AFDabsred<-AFDabs[!(AFDabs$AF_Plot=="D"|AFDabs$Metrics_Plot=="C"),]
NIELSENred<-NIELSEN[!(NIELSEN$AF_Plot=="D"|NIELSEN$Metrics_Plot=="C"),]
FLKred<-FLK[!(FLK$AF_Plot=="D"|FLK$Metrics_Plot=="C"),]
DDred<-DD[!(DD$AF_Plot=="D"|DD$Metrics_Plot=="C"),]

test<-Reduce(function(x,y) merge(x,y,by="Gene"),list(FSTred,DXYred,AFDabsred,NIELSENred,FLKred,DDred))

MZGenes<-unique(test$Gene)


FSTKK<-read.xlsx2("Genes_01percent_KletKowa.xlsx",1)
DXYKK<-read.xlsx2("Genes_01percent_KletKowa.xlsx",2)
AFDabsKK<-read.xlsx2("Genes_01percent_KletKowa.xlsx",3)
NIELSENKK<-read.xlsx2("Genes_01percent_KletKowa.xlsx",4)
VARLDKK<-read.xlsx2("Genes_01percent_KletKowa.xlsx",5)
FLKKK<-read.xlsx2("Genes_01percent_KletKowa.xlsx",6)

FSTKKred<-FSTKK[!(FSTKK$AF_Plot=="D"|FSTKK$Metrics_Plot=="C"),] 
DXYKKred<-DXYKK[!(DXYKK$AF_Plot=="D"|DXYKK$Metrics_Plot=="C"),]
AFDabsKKred<-AFDabsKK[!(AFDabsKK$AF_Plot=="D"|AFDabsKK$Metrics_Plot=="C"),]
NIELSENKKred<-NIELSENKK[!(NIELSENKK$AF_Plot=="D"|NIELSENKK$Metrics_Plot=="C"),]
FLKKKred<-FLKKK[!(FLKKK$AF_Plot=="D"|FLKKK$Metrics_Plot=="C"),]

testKK<-Reduce(function(x,y) merge(x,y,by="Gene"),list(FSTKKred,DXYKKred,AFDabsKKred,NIELSENKKred,FLKKKred))

KKGenes<-unique(testKK$Gene)

All1percentwoVarLD<-rbind(Thalgenes_described_Perc1Fst_lyratagenes_df_w_orthogroup,Thalgenes_described_Perc1DD_lyratagenes_df_w_orthogroup,Thalgenes_described_Perc1Nielsen_lyratagenes_df_w_orthogroup,Thalgenes_described_Perc1Dxy_lyratagenes_df_w_orthogroup,Thalgenes_described_Perc1Flk_lyratagenes_df_w_orthogroup,Thalgenes_described_Perc1AFDabs_lyratagenes_df_w_orthogroup)
write.table(unique(All1percentwoVarLD$Gene),"KletKowa_1percentoutliers_allwoVarLD.csv", sep=";", row.names=F,quote=F)




