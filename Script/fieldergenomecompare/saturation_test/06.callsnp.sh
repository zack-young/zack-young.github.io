#!/usr/bin/env bash
# Guo, Weilong; guoweilong@126.com; 2017-10-24
#for num in  3 5 6 7 8 9 11 13 ;do  #05_sample 10_sample;do # 25_sample 50_sample 75_sample 85_sample;do
#echo $CHR
num=9
chrlst=chr1A.1,chr1A.2,chr1B.1,chr1B.2,chr1D.1,chr1D.2,chr2A.1,chr2A.2,chr2B.1,chr2B.2,chr2D.1,chr2D.2,chr3A.1,chr3A.2,chr3B.1,chr3B.2,chr3D.1,chr3D.2,chr4A.1,chr4A.2,chr4B.1,chr4B.2,chr4D.1,chr4D.2,chr5A.1,chr5A.2,chr5B.1,chr5B.2,chr5D.1,chr5D.2,chr6A.1,chr6A.2,chr6B.1,chr6B.2,chr6D.1,chr6D.2,chr7A.1,chr7A.2,chr7B.1,chr7B.2,chr7D.1,chr7D.2
CHR=`echo ${chrlst} | tr "," " "`
DIR=`echo sample_depth_$num`
#echo $DIR
parallel -j procfile --results ${DIR} sh csp.sh ::: $(eval echo ${CHR}) ::: ${DIR}
#done
