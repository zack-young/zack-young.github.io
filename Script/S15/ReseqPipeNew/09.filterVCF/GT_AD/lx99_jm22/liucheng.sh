#!/bin/bash
set -x
for num in {1..7};do
    for i in {A,B,D}; do
       #cat ../MASK/chr${num}A_combine_total chr${num}A.lx99_density_graphy_gene|sort -nk2> chr${num}A_graphy
        bcftools view /data/user/yangzz/mapping/08.mergeGVCF/lx99_jm22/chr${num}${i}.ann.bcf.gz| bedtools intersect -a stdin -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/combine_snp_high -sorted -header |bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t[%GT\t%AD\t%GQ\t]\n' | /data/user/yangzz/mapping/09.filterVCF/GT_AD/muti_func_snp_compare.py --chromosome /data/user/yangzz/mapping/09.filterVCF/GT_AD/MASK/lx99_jm22_mask/chr${num}${i}.allCNV_use --two_diff on  --DP_low 3 --DP_high 20 --GQ_sample 8  -1 6 -2 9 -b 1000000 | awk '!a[$0]++{print}' > /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/annv1p1/chr${num}${i}.ann.high &
        bcftools view /data/user/yangzz/mapping/08.mergeGVCF/lx99_jm22/chr${num}${i}.ann.bcf.gz| bedtools intersect -a stdin -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/combine_snp_mid -sorted -header |bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t[%GT\t%AD\t%GQ\t]\n' | /data/user/yangzz/mapping/09.filterVCF/GT_AD/muti_func_snp_compare.py --chromosome /data/user/yangzz/mapping/09.filterVCF/GT_AD/MASK/lx99_jm22_mask/chr${num}${i}.allCNV_use --two_diff on  --DP_low 3 --DP_high 20 --GQ_sample 8  -1 6 -2 9 -b 1000000 | awk '!a[$0]++{print}' > /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/annv1p1/chr${num}${i}.ann.mid &

#        for a in {indels,snps}; do
#            bcftools view -v ${a} /data/user/yangzz/mapping/08.mergeGVCF/chr${num}${i}.ann_no_inter.bcf.gz| bedtools intersect -a stdin -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/combine_snp_high -sorted -header |bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t[%GT\t%AD\t%GQ\t]\n'> /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/annv1p1/snp_indel/chr${num}${i}_high.ann_${a} &
#            bcftools view -v ${a} /data/user/yangzz/mapping/08.mergeGVCF/chr${num}${i}.ann_no_inter.bcf.gz| bedtools intersect -a stdin -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/combine_snp_mid -sorted -header |bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t[%GT\t%AD\t%GQ\t]\n'> /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/annv1p1/snp_indel/chr${num}${i}_mid.ann_${a} &
    wait_all

    done
done
wait

for num in {1..7};do
    for i in {A,B,D}; do
        bcftools view /data/user/yangzz/mapping/08.mergeGVCF/lx99_jm22/chr${num}${i}.ann.bcf.gz|bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t[%GT\t%AD\t%GQ\t]\n' |/data/user/yangzz/mapping/09.filterVCF/GT_AD/muti_func_snp_compare.py --chromosome /data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/NZ10_YM05/chr${num}${i}.NZ10_YM05unmatch_homo_snp_level --two_diff on  --DP_low 3 --DP_high 20 --GQ_sample 8  -1 6 -2 9 -b 1000000 | awk '!a[$0]++{print}' > /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/annv1p1/chr${num}${i}.ann.all &
#        wait_all
#    awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a == 1) print a[$2]"\t"$4}' ../jm22_snp_indel_analysis/snp/chr${num}${i}_snp_altLEVEL ../lx99_snp_indel_analysis/snp/chr${num}${i}_snp_altLEVEL|sed '1d' > region/chr${num}${i}_solo_level_1.1
    #awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a == 1) print a[$2]"\t"$4}' region/chr${num}${i}_solo_level_1.1 region/chr${num}${i}_snp_level  > region/chr${num}${i}_combine_solo_plus_1.1
#     sed -i "1iCHR	start	end	alt_ratio_level_jm22	alt_ratio_level_lx99" region/chr${num}${i}_solo_level_1.1 
#     sed -i "1iCHR	start	end	alt_ratio_level_jm22	alt_ratio_level_lx99	diff_homosnp_ratio_level" region/chr${num}${i}_combine_solo_plus_1.1

    #../filtlevel.py --compare on -i region/chr${num}${i}_combine_solo_plus_1.1 >  region/chr${num}${i}_count_solo_plus

    done
done
#cat region/chr*_count_solo_plus > combine_count_solo_plus


#for num in {1..7};do
#    grep -F -f lx99_mask/chr${num}A.mask_CNV_test jm22_mask/chr${num}A.mask_CNV_test | sort -n -k3 | uniq > chr${num}A.mask_filtered_combine_test
#    grep -F -f lx99_mask/chr${num}B.mask_CNV_test jm22_mask/chr${num}B.mask_CNV_test | sort -n -k3 | uniq > chr${num}B.mask_filtered_combine_test
#    grep -F -f lx99_mask/chr${num}D.mask_CNV_test jm22_mask/chr${num}D.mask_CNV_test | sort -n -k3 | uniq > chr${num}D.mask_filtered_combine_test
#done
