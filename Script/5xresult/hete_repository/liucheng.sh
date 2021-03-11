#!/usr/bin/sh
for i in {1..7};do
    for a in {A,B,D};do
        echo "${i}${a}"
        ./muti_snp_hete_parent.py --DP_density on -i ~/mapping/09.filterVCF/GT_AD/lx99_jm22/chr${i}${a}_snp.gt -1 5 --chromosome ~/mapping/09.filterVCF/GT_AD/tmp/NZ10/chr${i}${a}.mask_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 > chr${i}${a}.DP_ratio &
        wait_all
    done
done
