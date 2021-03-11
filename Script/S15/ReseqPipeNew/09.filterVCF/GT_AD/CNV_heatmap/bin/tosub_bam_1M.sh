#!/usr/bin/env bash
set -euxo pipefail

while read i;do
  if [[ ! -d "/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/${i}" ]]; then
    mkdir /data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/${i}
  fi
  while read j;do
    sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/script.sh $i $j &
    wait_all
  done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/chrlst
  #wait_all
done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/rest_cultivar
# norm
#while read i;do
#  sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/script_1M.sh $i  &
#  wait_all
#done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt
