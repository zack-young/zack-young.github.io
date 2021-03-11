#!/usr/bin/env bash
set -euxo pipefail

while read i;do
  while read j;do
    sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/script.sh $i $j &
    wait_all
  done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/chrlst
done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt
# norm
while read i;do
  sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/normDPbySample.sh $i "1k.DP" "norm" 4 &
  wait_all
done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt
