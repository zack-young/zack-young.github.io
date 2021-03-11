#!/usr/bin/bash

SM=$1
path_SM="/data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/${1}_${5}"
OldSuffix=$2
NewSuffix=$3
DPcol=$4

# normlize DP coverage by sample
sum=$(gawk -vc=$DPcol '{sum+=$c} END{print sum}' ${path_SM}/*.${OldSuffix})
line=$(wc -l ${path_SM}/*.${OldSuffix}|tail -n 1|gawk '{print $1}')
ave=$(perl -e "print $sum/$line")

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  gawk -vOFS="\t" -vave=$ave -vc=$DPcol '{print $c/ave}' ${path_SM}/${CHR}.1.${OldSuffix} > ${path_SM}/${CHR}.${NewSuffix}
  gawk -vOFS="\t" -vave=$ave -vc=$DPcol '{print $c/ave}' ${path_SM}/${CHR}.2.${OldSuffix} >> ${path_SM}/${CHR}.${NewSuffix}
  sed -i '1i '$SM'' ${path_SM}/${CHR}.${NewSuffix}
done
