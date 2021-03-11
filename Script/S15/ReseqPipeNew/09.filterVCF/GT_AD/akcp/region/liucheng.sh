#!/usr/bin/env bash
for num in {1..7};do
    for i in {A,B,D}; do
        bcftools view /data/user/yangzz/mapping/08.mergeGVCF/chr${num}${i}.ann.bcf.gz|bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t[%GT\t%AD\t%GQ\t]\n' > chr${num}${i}.ann.gt &
#/data/user/yangzz/mapping/09.filterVCF/GT_AD/muti_func_snp_compare.py --two_diff on  --DP_low 3 --DP_high 20 --GQ_sample 8  -1 6 -2 9 -b 1000000 > /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/annv1p1/chr${num}${i}_${a}.ann.all & 
    done
done
