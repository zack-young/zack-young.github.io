#!/usr/bin/env bash
# Guo, Weilong; guoweilong@126.com; 2017-10-24

for BAM in ../05.mergeAsplit/$1/*.dedup.bam; do
  CHR=`basename $BAM | sed s/.dedup.bam//g`
  echo $CHR
  #wait_yhq
  sh csp.sh $CHR $1 2>&1 &
  wait_all
done
