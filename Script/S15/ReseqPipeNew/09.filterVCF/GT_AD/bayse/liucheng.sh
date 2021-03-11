#!/bin/bash

for num in {1..7}; do
   path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
  ./source_define.py -b 1000000 -i ${path}T123_combine_AD_mask/chr${num}A.gt_dp -1 5 -2 9 -s 7 -n ${num}A -o chr${num}A.bayes &
  ./source_define.py -b 1000000 -i ${path}T123_combine_AD_mask/chr${num}B.gt_dp -1 5 -2 9 -s 7 -n ${num}B -o chr${num}B.bayes &
  ./source_define.py -b 1000000 -i ${path}T123_combine_AD_mask/chr${num}D.gt_dp -1 5 -2 9 -s 7 -n ${num}D -o chr${num}D.bayes &
done
