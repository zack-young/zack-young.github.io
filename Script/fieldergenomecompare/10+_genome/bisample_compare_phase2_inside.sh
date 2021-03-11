#!/bin/bash
lines=$1
SAMPLE1=`echo $lines|awk -F ',' '{print $1}'`
SAMPLE2=`echo $lines|awk -F ',' '{print $2}'`
for num in {1..7}; do
  for i in {A,B,D}; do
      a=`du -h ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV`
      b=`echo $a |awk '{print $1}'`
      if [ "${b}" == 0 ];then
        awk '{print $0}'   ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.density  > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.delcnv_density
      else
        awk 'NR==FNR{a[$2]=$2;next}NR>FNR{if($2 in a ==0) print $0}'  ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV  ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.density  > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.delcnv_density
      fi

  done
done
