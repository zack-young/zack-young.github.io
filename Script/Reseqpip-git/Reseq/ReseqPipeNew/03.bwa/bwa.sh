#!/usr/bin/env bash
set -euxo pipefail

# Guo, Weilong; guoweilong@126.com; 2017-10-19

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start BWA"


DIR=`dirname $PWD`
SM=`basename $DIR`
#
ID=$1
R1="/data/user/yangzz/mapping/Reseqpip-git/Reseq/ReseqPipeNew/02.trimmomatic/$2/${ID}_R1.clean.fq.gz"
R2="/data/user/yangzz/mapping/Reseqpip-git/Reseq/ReseqPipeNew/02.trimmomatic/$2/${ID}_R2.clean.fq.gz"
#
bwa mem -t 20 \
  -R '@RG\tID:'${SM}'\tLB:'${SM}'\tPL:ILLUMINA\tSM:'${SM} \
  /data/index/wheat/IWGSCv2/bwa/iwgsc_refseqv2.0_all_chromosomes_parts.fa \
  ${R1} ${R2} | \
  samtools view -Sbh - > ${ID}.bam
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] To sort bam file"
#
samtools sort -@ 20 -m 48G -o ${ID}.sort.bam ${ID}.bam
#
samtools index ${ID}.sort.bam
#
samtools flagstat ${ID}.sort.bam > ${ID}.flagstat
#
rm ${ID}.bam*
#
wait
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish BWA"
#

