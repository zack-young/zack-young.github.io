#!/usr/bin/env bash
#set -euxo pipefail

line=$1
CHR=$2
time_date=$3
# calculate DP coverage for each sample & chr
bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${CHR}.1k.bed -b /data2/user2/zhangll/workdata/public/Elytrigia_repens/${CHR}.*.bam -counts > /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/${line}_${time_date}/${CHR}.1k.DP
