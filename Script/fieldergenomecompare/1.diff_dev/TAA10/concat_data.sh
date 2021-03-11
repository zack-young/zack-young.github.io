#!/bin/bash

#if false;then
#for i in {1..595};do
#  datamash transpose < chr1A.${i}.cluslabel
#done >> chr1A.allsample.cluslabel

#((echo "CHROM\nstart\nend"; cat ../WGS_samplelist.txt)| datamash transpose; ( paste chr1A.1M.bed chr1A.allsample.cluslabel) ) > chr1A.allsample.cluslabel.fomt
#fi

#Rscript /data2/rawdata2/tetraintro/SNPClusterInBin/Clustmap.R -k 10 -c metadata.txt -m namemap.txt -i chr1A.allsample.cluslabel.fomt -H 10 -W 11 -t chr1A -w 1 -o chr1A.pdf
chr=$1
path=$2
suffix=$3
BED=$4
num=`cat ${BED}|grep $chr | wc -l`
for i in `seq  1 $num`;do
echo "${path}/${chr}.${i}.${suffix}"
done|tr '\n' ' '|sed s/' '$//g
