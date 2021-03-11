#!/usr/bin/env bash
# Guo, Weilong; guoweilong@126.com; 2017-10-24
for DIR in 05_sample 10_sample 25_sample 50_sample 75_sample 85_sample;do
for BAM in ${DIR}/*.dedup.bam; do
  CHR=`basename $BAM | sed s/.dedup.bam//g`
  #echo $CHR
  sh csp.sh $CHR $DIR > ${DIR}/${CHR}.o 2>${DIR}/${CHR}.e &
wait_all
done
done
