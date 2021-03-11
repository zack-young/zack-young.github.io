#!/usr/bin/env bash
set -xuo pipefail
# Guo, Weilong; guoweilong@126.com; 2017-10-19
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start Split"


#for ID in `for F in ../01*/R1_??.gz; do basename $F; done | cut -c4-5 | tr '\n' ' '`; do
for ID in `while read line; do echo $line|tr -d '\n'; done < $1`; do 
  echo $ID
  R1="R1_"${ID}.gz
  R2="R2_"${ID}.gz
  set +e
  mkdir ${2}
  set +e
  cd ${2}
  (java \
     -jar /home/software/trinityrnaseq-Trinity-v2.4.0/trinity-plugins/Trimmomatic-0.36/trimmomatic-0.36.jar \
     PE -phred33 -threads 1 \
     ../../01.split/${2}/${R1} ../../01.split/${2}/${R2} \
     ${ID}_R1.clean.fq.gz ${ID}_R1.unpaired.fq.gz ${ID}_R2.clean.fq.gz ${ID}_R2.unpaired.fq.gz \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;
   zcat ../../01.split/${2}/${R1} | gawk '{T++;} END{print T/4;}' > ../../01.split/${2}/${ID}.count ;
   rm ../../01.split/${2}/${R1} ../../01.split/${2}/${R2} ;
   zcat ${ID}_R1.clean.fq.gz | gawk '{T++;} END{print T/4;}' > ${ID}.count )&
  cd ..
done

wait

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish Split"

#cd ../03.*/
#sh 03.*.sh

