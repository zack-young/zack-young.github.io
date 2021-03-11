#!/bin/bash

#for item in {1..3};do
#    a=`expr $item / 2 - 1`
#    for num in {1..7};do
#       cat T${item}_mask/chr${num}A.mask_CNV | awk '{print $1"\t"$3"\t"$4}' > T${item}_mask/chr${num}A.mask_CNVfiltered
#       cat T${item}_mask/chr${num}B.mask_CNV | awk '{print $1"\t"$3"\t"$4}' > T${item}_mask/chr${num}B.mask_CNVfiltered
#       cat T${item}_mask/chr${num}D.mask_CNV | awk '{print $1"\t"$3"\t"$4}' > T${item}_mask/chr${num}D.mask_CNVfiltered
##       awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""'${a}'""\t"$6"\t"$7"\t"$8}' s${a}_AD/chr${num}D.graphy > s${a}_AD/chr${num}D.facet_ne
#    done
#    wait
#done
for num in {1..7};do
    for i in {A,B,D};do
        #cat lx99_mask/chr${num}A.mask_CNV_test jm22_mask/chr${num}A.mask_CNV_test | sort -nk3 -u -n|cut -f1,3,4,5> chr${num}A.mask_filtered_combine_DP
        #cat chr${num}A.mask_filtered_combine_test /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/CNV/combine_filter_bin |grep "chr${num}A" | sort -nk3 -u > chr${num}A.mask_filtered_combine_thre
        #grep 'CNV' chr${num}A.mask_filtered_combine_thre > chr${num}A.CNV
        #awk '{print $1"\t"$2"\t"$3"\t""both"}' chr${num}${i}.bothCNV > chr${num}${i}.bothCNV_use
        #awk '{print $1"\t"$2"\t"$3"\t""jm22"}' lx99_jm22_mask/chr${num}${i}_jm22have > lx99_jm22_mask/chr${num}${i}_jm22have_use
        #awk '{print $1"\t"$2"\t"$3"\t""lx99"}' lx99_jm22_mask/chr${num}${i}_lx99have > lx99_jm22_mask/chr${num}${i}_lx99have_use
        cat chr${num}${i}.bothCNV_use lx99_jm22_mask/chr${num}${i}_jm22have_use lx99_jm22_mask/chr${num}${i}_lx99have_use|sort -nk2 > chr${num}${i}.allCNV_use
        #cat jm22_mask/chr${num}${i}.mask_CNV lx99_mask/chr${num}${i}.mask_CNV |sort -nk2 -u > chr${num}${i}_all_CNV
        #grep -F -f /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/MASK/lx99_mask/chr${num}${i}.mask_CNV /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/MASK/jm22_mask/chr${num}${i}.mask_CNV | sort -nk2 -u> chr${num}${i}.bothCNV
        #awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a == 0) print $0}' jm22_mask/chr${num}${i}.mask_CNV lx99_mask/chr${num}${i}.mask_CNV > lx99_jm22_mask/chr${num}${i}_lx99have
        #awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a == 0) print $0}' lx99_mask/chr${num}${i}.mask_CNV jm22_mask/chr${num}${i}.mask_CNV > lx99_jm22_mask/chr${num}${i}_jm22have
    done
done
#cat lx99_jm22_mask/chr*_lx99have > lx99_jm22_mask/combine_lx99have
#cat lx99_jm22_mask/chr*_jm22have > lx99_jm22_mask/combine_jm22have

#for num in {1..7};do
#    grep -F -f lx99_mask/chr${num}A.mask_CNV_test jm22_mask/chr${num}A.mask_CNV_test | sort -n -k3 -u > chr${num}A.mask_filtered_combine_test
#done
