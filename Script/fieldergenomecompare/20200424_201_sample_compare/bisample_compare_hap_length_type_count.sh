#!/bin/bash
#set -euxo pipefail
CHR=$1
line=$2
SAMPLE1=`echo $line|awk -F ',' '{print $1}'`
SAMPLE2=`echo $line|awk -F ',' '{print $2}'`

if true;then
num=`sort -nk2,2 /data3/user3/wangwx/projs/HMM_for_yzz_comp/201212/201_sample_compare_HMMdata/${SAMPLE1}_${SAMPLE2}/${CHR}.homo_undefined_snp_level.sorted.v4|grep -E "low"|wc -l` 
if [  $num != 0 ];then
sort -nk2,2 /data3/user3/wangwx/projs/HMM_for_yzz_comp/201212/201_sample_compare_HMMdata/${SAMPLE1}_${SAMPLE2}/${CHR}.homo_undefined_snp_level.sorted.v4|grep -E "low"|bedtools merge -i - -c 4 -o distinct > ~/mapping/fieldergenomecompare/20200424_201_sample_compare/${SAMPLE1}_${SAMPLE2}/${CHR}.hmm.v4.low_both_merge.txt
else
rm -f ~/mapping/fieldergenomecompare/20200424_201_sample_compare/${SAMPLE1}_${SAMPLE2}/${CHR}.hmm.v4.low_both_merge.txt
touch ~/mapping/fieldergenomecompare/20200424_201_sample_compare/${SAMPLE1}_${SAMPLE2}/${CHR}.hmm.v4.low_both_merge.txt
fi
fi

