#!/usr/bin/env bash
set -xeuo pipefail

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start Split"

batch=~/mapping/09.filterVCF/GT_AD/lx99_jm22/fastq
for i in ${batch}/*_1*;do
    SM=`basename $i|cut -d '_' -f1`
    set +e
    mkdir $SM;
    set -e
    cd $SM;
    zcat $i | split -l 100000000 - R1_ &
    zcat ${batch}/${SM}*_2* | split -l 100000000 - R2_ &
    wait
    for F in R?_??; do
        gzip $F &
        wait_all
    done
    wait
    echo`ls R1_* | sed s/R1_//g | sed s/.gz//g| tr '\n' ' '` > ID.list
    cd ..
    #GZrsync ${SM} &
done
#
wait
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish Split"
cd ../02.*/
nohup sh trimmomatic.sh ../01.split/C4-2/ID.list C4-2 > C4-2.out  2>&1 &
nohup sh trimmomatic.sh ../01.split/LX99-1-1/ID.list LX99-1-1 > LX99-1-1.out  2>&1 &
