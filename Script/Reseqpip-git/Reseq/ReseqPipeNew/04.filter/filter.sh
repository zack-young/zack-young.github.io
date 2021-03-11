#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19
#
#source /WORK/app/osenv/ln1/set2.sh
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start filter"
#
WP="/WORK/pp192/"
#
set +e
mkdir $1
set -e
for BAM in ../03.*/$1/??.sort.bam; do
  echo $BAM
  ID=`basename ${BAM} | sed s/.sort.bam//g`
  (bamtools filter \
    -in ${BAM}  \
    -out $1/${ID}.filter.bam \
    -forceCompression -script filter.list;
   samtools index $1/${ID}.filter.bam > $1/${ID}.filter.bam.bai;
   samtools flagstat $1/${ID}.filter.bam > $1/${ID}.flagstat
   )&
#   rm ${BAM})&
done
#
wait
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish Filter"
#
#cd ../05*/
#sh ./05.mergeAsplit.sh
#
