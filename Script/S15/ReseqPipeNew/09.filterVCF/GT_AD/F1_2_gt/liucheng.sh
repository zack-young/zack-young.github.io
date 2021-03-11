#!/bin/bash
#for num in {1..7};do
#    awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' chr${num}A.2_f1_2.gt > chr${num}A.2_f1_2.gt.log
#    awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' chr${num}B.2_f1_2.gt > chr${num}B.2_f1_2.gt.log
#    awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' chr${num}D.2_f1_2.gt > chr${num}D.2_f1_2.gt.log
#    cat chr${num}A.1_f1_2.gt chr${num}A.2_f1_2.gt.log > chr${num}A.gt 
#    cat chr${num}B.1_f1_2.gt chr${num}B.2_f1_2.gt.log > chr${num}B.gt 
#    cat chr${num}D.1_f1_2.gt chr${num}D.2_f1_2.gt.log > chr${num}D.gt 
#done
#rm *gt.log
#path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
#for num in {1..7};do
#    ./snp_density.py -b 1000000 -i ${path}F1_AD/chr${num}A.gt -1 5 > ${path}F1_2_gt/chr${num}A.f1_snp_count &
#    ./snp_density.py -b 1000000 -i ${path}F1_AD/chr${num}B.gt -1 5 > ${path}F1_2_gt/chr${num}B.f1_snp_count &
#    ./snp_density.py -b 1000000 -i ${path}F1_AD/chr${num}D.gt -1 5 > ${path}F1_2_gt/chr${num}D.f1_snp_count &
#    ./snp_density.py -b 1000000 -i ${path}F2_AD/chr${num}A.gt -1 5 > ${path}F1_2_gt/chr${num}A.f2_snp_count &
#    ./snp_density.py -b 1000000 -i ${path}F2_AD/chr${num}B.gt -1 5 > ${path}F1_2_gt/chr${num}B.f2_snp_count &
#    ./snp_density.py -b 1000000 -i ${path}F2_AD/chr${num}D.gt -1 5 > ${path}F1_2_gt/chr${num}D.f2_snp_count &
#done

path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
for num in {1..7};do
    ./ratio_filter.py -b 1000000 -i ${path}F1_2_gt/chr${num}A.f1_threshold_DP_GQ -1 5 > ${path}F1_2_gt/chr${num}A.f1_hete_count &
    ./ratio_filter.py -b 1000000 -i ${path}F1_2_gt/chr${num}B.f1_threshold_DP_GQ -1 5 > ${path}F1_2_gt/chr${num}B.f1_hete_count &
    ./ratio_filter.py -b 1000000 -i ${path}F1_2_gt/chr${num}D.f1_threshold_DP_GQ -1 5 > ${path}F1_2_gt/chr${num}D.f1_hete_count &
    ./ratio_filter.py -b 1000000 -i ${path}F1_2_gt/chr${num}A.f2_threshold_DP_GQ -1 5 > ${path}F1_2_gt/chr${num}A.f2_hete_count &
    ./ratio_filter.py -b 1000000 -i ${path}F1_2_gt/chr${num}B.f2_threshold_DP_GQ -1 5 > ${path}F1_2_gt/chr${num}B.f2_hete_count &
    ./ratio_filter.py -b 1000000 -i ${path}F1_2_gt/chr${num}D.f2_threshold_DP_GQ -1 5 > ${path}F1_2_gt/chr${num}D.f2_hete_count &
done


#for num in {1..7}; do     # count different snp between lx99 and jm22
#  path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
#  ./DP_GQ_filter.py -b 1000000 -i chr${num}A.gt -1 5 -o chr${num}A.filtered_DP_GQ &
#  ./DP_GQ_filter.py -b 1000000 -i chr${num}B.gt -1 5 -o chr${num}B.filtered_DP_GQ &
#  ./DP_GQ_filter.py -b 1000000 -i chr${num}D.gt -1 5 -o chr${num}D.filtered_DP_GQ &
#done

#for num in {1..7}; do     # count different snp between lx99 and jm22
#  path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
#  ./threshold-filter.py -b 1000000 -i chr${num}A.filtered_DP_GQ -1 5 -o chr${num}A.f1_threshold_DP_GQ &
#  ./threshold-filter.py -b 1000000 -i chr${num}B.filtered_DP_GQ -1 5 -o chr${num}B.f1_threshold_DP_GQ &
#  ./threshold-filter.py -b 1000000 -i chr${num}D.filtered_DP_GQ -1 5 -o chr${num}D.f1_threshold_DP_GQ &
#  ./threshold-filter.py -b 1000000 -i chr${num}A.filtered_DP_GQ -1 8 -o chr${num}A.f2_threshold_DP_GQ &
#  ./threshold-filter.py -b 1000000 -i chr${num}B.filtered_DP_GQ -1 8 -o chr${num}B.f2_threshold_DP_GQ &
#  ./threshold-filter.py -b 1000000 -i chr${num}D.filtered_DP_GQ -1 8 -o chr${num}D.f2_threshold_DP_GQ &
#done

