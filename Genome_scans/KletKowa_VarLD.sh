#!/bin/bash

#run from pflaphy-gscan/GS_KletKowa_lyrata

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N VarLD_KletKowa
#$ -l vf=2G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

cut -f1-13 VarLD_KletKowa_arenosa.csv > VarLD_Klet.table
cut -f1-3,14-20 VarLD_KletKowa_arenosa.csv > VarLD_Kowa.table

chr=1
grep "scaffold_1[^[:digit:]]"  VarLD_Klet.table | cut -f2-13 > File1.txt
grep "scaffold_1[^[:digit:]]"  VarLD_Kowa.table | cut -f2-10 > File2.txt

java -XX:ParallelGCThreads=2 -Xmx50g -jar ../software/rgenetics-1.0.jar -p VarLD -n 25 -f 0 -m 0 -o KletKowa_VarLD.out File1.txt File2.txt 

cp KletKowa_VarLD.out KletKowa_VarLD.out.table
length=$(wc -l KletKowa_VarLD.out | cut -d' ' -f1)
printf "scaffold_"$chr"\n%.0s" $(seq 1 $length) > Chr.txt

for chr in $(seq 2 8)
do
	grep "scaffold_"$chr"[^[:digit:]]" VarLD_Klet.table | cut -f2-13 > File1.txt
	grep "scaffold_"$chr"[^[:digit:]]" VarLD_Kowa.table | cut -f2-10 > File2.txt
	
	File1length=$(wc -l File1.txt | cut -f1 -d ' ')
        File2length=$(wc -l File2.txt | cut -f1 -d ' ')
        echo "Chr"$chr"_File1_"$File1length"_File2_"$File2length >> File12length.txt
        if [ $File1length -ne "0" ]
        then
                java -XX:ParallelGCThreads=2 -Xmx50g -jar ../software/rgenetics-1.0.jar -p VarLD -n 25 -f 0 -m 0 -o KletKowa_VarLD.out File1.txt File2.txt
        	tail -n +2 KletKowa_VarLD.out > KletKowa_VarLD.txt
        	cat KletKowa_VarLD.txt >> KletKowa_VarLD.out.table
        	length=$(wc -l KletKowa_VarLD.txt | cut -d' ' -f1)
        	printf "Chr_"$chr"Length_"$length"\n" >> Length.txt
        	if [ $length -ne 0 ]
        	then
                	printf "scaffold_"$chr"\n%.0s" $(seq 1 $length) >> Chr.txt
                	printf $chr" " >> Test.txt
        	fi
	fi
	> KletKowa_VarLD.out
done
#run from VarLD folder

