#take numbers from Venn diagram

#phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
#where q=size of overlap-1; m=number of upregulated genes in experiment #1; 
#n=(total number of genes on platform-m); k=number of upregulated genes in experiment #2.


#‘pick m balls from a jar of n balls, then repeat by picking k balls: what is the significance of exactly q overlap. 

#halleri
q=35-1
m=840+32+7+2+9+3
k=204+32+3+3+2
n=33221-(840+32+7+2+9+3)

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#7.313168e-16

sum(dhyper((q+1):k,m,n-m,k))
fisher.test(matrix(c(n-m-k+q+1,k-q+1,m-q+1,q+1), nrow=2), alternative="greater")
#        Fisher's Exact Test for Count Data

#data:  matrix(c(n - m - k + q + 1, k - q + 1, m - q + 1, q + 1), nrow = 2)
#p-value = 2.249e-15
#alternative hypothesis: true odds ratio is greater than 1
#95 percent confidence interval:
# 4.325351      Inf
#sample estimates:
#odds ratio 
#  6.022049 
#phyper(q+1, m, n, k, lower.tail = F, log.p = FALSE)
#[1] 1.127444e-16




#arenosa
q=6-1
m=146+4+2+9+3+2
k=165+4+7+2+3
n=33221-(146+4+2+9+3+2)
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.0003107269

#MZ halleri arenosa
q=3-1
m=204+32+3+3+2
k=165+4+7+2+3
n=33221-(204+32+3+3+2)
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.1486956

#KK halleri arenosa
q=9+3+2-1
m=840+32+7+2+9+3
k=146+4+2+9+3+2
n=33221-(840+32+7+2+9+3)
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.0001686143

#KK halleri MZ arenosa
q=7+2-1
m=840+32+7+2+9+3
k=165+4+7+2+3
n=33221-(840+32+7+2+9+3)
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.05659832

#MZ halleri KK arenosa 
q=3+2-1
m=204+32+3+3+2
k=146+4+2+9+3+2
n=33221-(204+32+3+3+2)
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.007796879

genelist<-read.table("Alyrata_107_v1.0.annotation.v2_Feb2015.gff.txt",sep="\t",quote="")
#V4=start_pos
#V5=end_pos
names(genelist)<-c("Scaffold","Version","Type","Start_pos","End_pos","Score","Strand","Frame","ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)
IDlist1<-strsplit(as.character(genelist[,9]),";")
IDs1<-matrix(unlist(IDlist1),ncol=3,byrow=TRUE)
IDlist2<-strsplit(as.character(IDs1[,1]),"=")
IDs2<-as.data.frame(matrix(unlist(IDlist2),ncol=2,byrow=TRUE))
genelist$ID<-as.character(IDs2[,2])

#33,221 genes

#GS only

#halleri
q=4+1-1
k=67+4+1+1
m=85+4+2+2+1
n=33221-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#2.101362e-06

#sum(dhyper((q+1):k,m,n-m,k))
fisher.test(matrix(c(n-m-k+q+1,k-q+1,m-q+1,q+1), nrow=2), alternative="greater")
#        Fisher's Exact Test for Count Data
#
#data:  matrix(c(n - m - k + q + 1, k - q + 1, m - q + 1, q + 1), nrow = 2)
#p-value = 2.703e-06
#alternative hypothesis: true odds ratio is greater than 1
#95 percent confidence interval:
# 9.697018      Inf
#sample estimates:
#odds ratio 
#  25.87131 

#arenosa
q=5-1
m=139+5+1+2
k=127+5+1+2
n=33221-m
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.0003461839

#MZ halleri arenosa
q=2-1
m=127+5+1+2
k=85+4+2+2+1
n=33221-m
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.05622982

#KK halleri arenosa
q=1-1
m=139+5+1+2
k=67+4+1+1
n=33221-m
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.276813

#KK halleri MZ arenosa
q=1-1
m=127+5+1+2
k=67+4+1+1
n=33221-m
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.2573849

#MZ halleri KK arenosa 
q=2+1-1
m=139+5+1+2
k=85+4+2+1+2
n=33221-m
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.008482792


#GS only and n set to 25,000

#halleri
q=4+1-1
k=67+4+1+1
m=85+4+2+2+1
n=25000-m

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#8.283806e-06

fisher.test(matrix(c(n-m-k+q+1,k-q+1,m-q+1,q+1), nrow=2), alternative="greater")
#        Fisher's Exact Test for Count Data

#data:  matrix(c(n - m - k + q + 1, k - q + 1, m - q + 1, q + 1), nrow = 2)
#p-value = 1.067e-05
#alternative hypothesis: true odds ratio is greater than 1
#95 percent confidence interval:
# 7.279524      Inf
#sample estimates:
#odds ratio 
#  19.41478 

#arenosa
q=5-1
m=139+5+1+2
k=127+5+1+2
n=25000-m
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.001233436

#MZ halleri arenosa
q=2-1
m=127+5+1+2
k=85+4+2+2+1
n=25000-m
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.09182121

#KK halleri arenosa
q=1-1
m=139+5+1+2
k=67+4+1+1
n=25000-m
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.3502255

#KK halleri MZ arenosa
q=1-1
m=127+5+1+2
k=67+4+1+1
n=25000-m
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.3268856

#MZ halleri KK arenosa 
q=2+1-1
m=139+5+1+2
k=85+4+2+1+2
n=25000-m
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.01809078


