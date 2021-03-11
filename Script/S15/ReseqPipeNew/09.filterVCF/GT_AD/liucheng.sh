#!/bin/bash
set -x
#for num in {1..7};do
#    for i in {A,B,D};do
        #./hete_define.py -b 1000000 -i S15_NZ10/chr${num}${i}_snp.gt --chromosome tmp/S15_2019-11-23224422/chr${num}${i}.mask_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 --hete_count on > hete_BanSheng/S15/chr${num}${i}_snp_hete_homo.count &
        #./hete_define.py -b 1000000 -i S60_NZ10/chr${num}${i}_snp.gt --chromosome tmp/S60_2019-11-23224505/chr${num}${i}.mask_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 --hete_count on > hete_BanSheng/S60/chr${num}${i}_snp_hete_homo.count &
        #./hete_define.py -b 1000000 -i S61_NZ10/chr${num}${i}_snp.gt --chromosome tmp/S61_2019-11-23224547/chr${num}${i}.mask_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 --hete_count on > hete_BanSheng/S61/chr${num}${i}_snp_hete_homo.count &
        #./hete_define.py -b 1000000 -i S62_NZ10/chr${num}${i}_snp.gt --chromosome tmp/S62_2019-11-23224741/chr${num}${i}.mask_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 --hete_count on > hete_BanSheng/S62/chr${num}${i}_snp_hete_homo.count &
        
#        ./hete_define.py -b 1000000 -i S15_NZ10/chr${num}${i}_snp.gt --chromosome tmp/S15_2019-11-23224422/chr${num}${i}.mask_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 --hete_define on > hete_BanSheng/S15/chr${num}${i}_snp_hete_homo.count &
#        ./hete_define.py -b 1000000 -i S60_NZ10/chr${num}${i}_snp.gt --chromosome tmp/S60_2019-11-23224505/chr${num}${i}.mask_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 --hete_define on > hete_BanSheng/S60/chr${num}${i}_snp_hete_homo.count &
#        ./hete_define.py -b 1000000 -i S61_NZ10/chr${num}${i}_snp.gt --chromosome tmp/S61_2019-11-23224547/chr${num}${i}.mask_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 --hete_define on > hete_BanSheng/S61/chr${num}${i}_snp_hete_homo.count &
#        ./hete_define.py -b 1000000 -i S62_NZ10/chr${num}${i}_snp.gt --chromosome tmp/S62_2019-11-23224741/chr${num}${i}.mask_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 --hete_define on > hete_BanSheng/S62/chr${num}${i}_snp_hete_homo.count &
#        wait_all
#    done
#done

for num in {1..7};do
    for i in {A,B,D};do
        bcftools view  -v snps -Ov ~/mapping/08.mergeGVCF/field_cultivar/chr${num}${i}.bcf.gz| bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%GQ\t]\n' > chr${num}${i}.snp_matrix &
        bcftools view  -v indels -Ov ~/mapping/08.mergeGVCF/field_cultivar/chr${num}${i}.bcf.gz| bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%GQ\t]\n' > chr${num}${i}.indel_matrix &
#        for item in {snp,indel};do
#        ./muti_func_snp_compare.py --density_graph_raw_vcf_one on -i XC01_XC02/chr${num}${i}_${item}.gt -1  8 --chromosome tmp/XC02_2019-12-11201602/chr${num}${i}.mask_CNV --DP_low 3 --DP_high 99 --GQ_sample 8 -b 1000000 --LEVEL "1.5,3.0" -o test/chr${num}${i}.XC02_${item}_density &
#        ./muti_func_snp_compare.py --density_graph_raw_vcf_one on -i XC01_XC02/chr${num}${i}_${item}.gt -1  5 --chromosome tmp/XC01_2019-12-11201602/chr${num}${i}.mask_CNV --DP_low 3 --DP_high 99 --GQ_sample 8 -b 1000000 --LEVEL "1.5,3.0" -o test/chr${num}${i}.XC01_${item}_density &
#    done
wait_all
done
done



#for num in {1..7};do
#    for i in {A,B,D};do
        #bcftools view -v snps ~/mapping/08.mergeGVCF/lx99_jm22/chr${num}${i}.ann.bcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]%ANN\n' - > lx99_jm22/chr${num}${i}_snp.gt &

#        ./muti_func_snp_filter.py -i lx99_jm22/chr${num}${i}_snp.gt -b 1000000 --chromsome tmp/NZ10_YM05/chr${num}${i}.NZ10_YM05unmatch_homo_snp_level --snp_indel_count_compare on --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 -2 8  -o lx99_jm22/snp/chr${num}${i}_ann_counts &
#        join -t $'\t' -1 2 -2 2 -o 1.1 -o 1.2 -o 1.3 -o 1.4 -o 2.4 tmp/CP06_CP11/chr${num}${i}.CP06_CP11unmatch_homo_snp_level   tmp/CP06_CP08/chr${num}${i}.CP06_CP08unmatch_homo_snp_level|sed '1c CHR	start	end	CP06_CP11	CP06_CP08'> CP06_CP08_CP11/chr${num}${i}CP06_CP11_CP08_combine_unmatch_level
        #./find_parents.py -i CP06_CP08_CP11/chr${num}${i}_snp.gt -b 1000000 --pre_define on --chromosome CP06_CP08_CP11/chr${num}${i}CP06_CP11_CP08_combine_unmatch_level --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 -2 11 -s 8 -o CP06_CP08_CP11/chr${num}${i}_parent_define &
        #./muti_func_snp_filter.py -i lx99_jm22/chr${num}${i}_indel.gt -b 1000000 --chromsome tmp/NZ10_YM05/chr${num}${i}.NZ10_YM05unmatch_homo_snp_level --snp_indel_count_compare on --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 -2 8  -o lx99_jm22/snp/chr${num}${i}_indel_counts 
        #grep -E 'mid_com|high_com' tmp/NZ10_YM05/chr${num}${i}.NZ10_YM05unmatch_homo_snp_level |awk 'BEGIN{sum=0}{sum+=$5}END{print sum}' > lx99_jm22/snp/chr${num}${i}_diff_snp_counts &
        #grep -E 'low_com' tmp/NZ10_YM05/chr${num}${i}.NZ10_YM05unmatch_homo_snp_level |awk 'BEGIN{sum=0}{sum+=$5}END{print sum}' > lx99_jm22/snp/chr${num}${i}_simi_snp_counts &
#    done
#done
#wait
#for i in {A,B,D};do
#    cat lx99_jm22/snp/chr*${i}_indel_counts > lx99_jm22/snp/combine_${i}indel_counts
    #cat lx99_jm22/snp/chr*_diff_snp_counts > lx99_jm22/snp/combine_diffsnp_counts
    #cat lx99_jm22/snp/chr*_simi_snp_counts > lx99_jm22/snp/combine_simisnp_counts
#done
#for item in {5..33..2};do
#for item in {5..9..2};do
#    a=`expr $item / 2 - 1`
#    for num in {1..7};do
#        ./filter.py -i chr${num}A_3_sample.gt -1 ${item} > T${a}_AD/chr${num}A.gt &  # filter unqualified snp (one of AD is 1 in heterogous ) 
#        ./filter.py -i chr${num}A.2_f1_2.gt -1 ${item} > F${a}_AD/chr${num}A.2.gt &
#        awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6}' F${a}_AD/chr${num}A.2.gt > F${a}_AD/chr${num}A.2.4b.gt &
#        awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6}' F${a}_AD/chr${num}B.2.gt > F${a}_AD/chr${num}B.2.4b.gt &
#        awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6}' F${a}_AD/chr${num}D.2.gt > F${a}_AD/chr${num}D.2.4b.gt &
#        cat F${a}_AD/chr${num}A.1.gt F${a}_AD/chr${num}A.2.4b.gt > F${a}_AD/chr${num}A.gt &
#        cat F${a}_AD/chr${num}B.1.gt F${a}_AD/chr${num}B.2.4b.gt > F${a}_AD/chr${num}B.gt &
#        cat F${a}_AD/chr${num}D.1.gt F${a}_AD/chr${num}D.2.4b.gt > F${a}_AD/chr${num}D.gt &

#       if [ ! -d "DP/T${a}_DP" ]; then
#            mkdir DP/T${a}_DP
#       fi
#        ./filter.py -i intial_gt/chr${num}D.2.gt -1 ${item} > s${a}_AD/chr${num}D.2.gt &
#        awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6}' s${a}_AD/chr${num}A.2.gt > s${a}_AD/chr${num}A.2.4b.gt &
#        cat s${a}_AD/chr${num}A.1.gt s${a}_AD/chr${num}A.2.4b.gt > s${a}_AD/chr${num}A.gt &
#       echo s${a}_AD/chr${num}A.graphy
#       ./graph.py -b 1000000 -n S${a} -1 5 -i s${a}_AD/chr${num}A.gt > s${a}_AD/chr${num}A.graphy &
#       echo s${a}_AD/chr${num}B.graphy
#       ./graph.py -b 1000000 -n S${a} -1 5 -i s${a}_AD/chr${num}B.gt > s${a}_AD/chr${num}B.graphy &
#       echo s${a}_AD/chr${num}D.graphy
#       ./graph.py -b 1000000 -n S${a} -1 5 -i s${a}_AD/chr${num}D.gt > s${a}_AD/chr${num}D.graphy &
#       ./unCNV_graph.py -b 1000000 -n S${a} -1 5 -i s${a}_AD/chr${num}A.gt > s${a}_AD/chr${num}A.graphy &
#       cat DP/s${a}_DP/*.DP > DP/s${a}_DP/combine_data 
#       mkdir graphy/s${a}_chr_graphy
#       awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""'${a}'""\t"$6"\t"$7"\t"$8}' s${a}_AD/chr${num}A.graphy > s${a}_AD/chr${num}A.facet_ne
#       find  -name "chr${num}A.graphy" -exec cat '{}' \;|sort -n -k 2 > graphy/chr${num}A.facet.combine_ne &
#./ratio_filter.py -b 1000000 -m 9.89 -l 9.23 -n S1 -1 5 -i s1_AD/chr1A.gt > filtered_ratio/s1_ratio/chr1A.ratio
#        wait_all
#    done
#    wait
#    cat DP/T${a}_DP/*.DP > DP/T${a}_DP/combine_data
#done
#wait
#cat DP/s1_DP/*.DP.new > DP/s1_DP/combine_data.new
#combine 0 combine_ne -0.1
#item=5
#a=123
#for num in {1..7};do
#    ./filter_multisample.py -i 3_sample_gt/chr${num}A_3_sample.gt -1 ${item} > T${a}_combine_AD_mask/chr${num}A.gt_dp &
#done

#for num in {1..7};do
#    ./lx99_jm22_filter_multisample.py -i chr${num}A_lx99_jm22.gt -1 5 > lx99_jm22/chr${num}A.gt_dp &
#done

#for num in {1..7}; do
#   path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
#  ./source_define.py -b 5000000 -i ${path}T123_combine_AD_mask/chr${num}B.gt -1 5 -2 7 -s 6 -n ${num}B -o ${path}T123_combine_AD_mask/benzi/chr${num}B.mask_filter &
#  ./source_define.py -b 5000000 -i ${path}T123_combine_AD_mask/chr${num}D.gt -1 5 -2 7 -s 6 -n ${num}D -o ${path}T123_combine_AD_mask/benzi/chr${num}D.mask_filter &
#done
#for a in {1..2}; do
#    for num in {1..7}; do
#        path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
#        ./density.py -b 1000000 -n F${a} -1 5 -i F${a}_AD/chr${num}A.gt > F${a}_AD/chr${num}A.density &
#    done
#done
#for num in {1..7}; do     # count different snp between lx987 and 3097
#  path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
  #CHR=`echo $VCF | sed s/.gt//g`
#  ./lx987_3097_different.py -b 1000000 -i ${path}T123_combine_AD_mask/chr${num}A.gt -1 5 -2 7 -s 6 -n ${num}A -o ${path}T123_combine_AD_mask/chr${num}A.lx987_3097 &
#done

#for num in {1..7}; do     # count different snp between lx99 and jm22
#  path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
# #CHR=`echo $VCF | sed s/.gt//g`
#  ./two_geno_diff.py -b 1000000 -i ${path}lx99_jm22/chr${num}A.gt_dp -1 5 -2 7  -n ${num}A -o ${path}lx99_jm22/chr${num}A.diff &
#done

#for num in {1..7}; do
#    if [ ! -d "MASK/lx99_mask" ]; then
#         mkdir MASK/lx99_mask
#    fi
#    if [ ! -d "MASK/jm22_mask" ]; then
#         mkdir MASK/jm22_mask
#    fi
#
#    sample1=/data2/rawdata2/readDepth/YM05
#    sample2=/data2/rawdata2/readDepth/NZ10
#    for i in {A,B,D};do
#        ./muti_func_snp_compare.py -b 1000000 -i /data2/rawdata2/readDepth/NZ10/chr${num}${i}.1M.norm --mask_cnv on -o /data/user/yangzz/mapping/09.filterVCF/GT_AD/MASK/jm22_mask/chr${num}${i}.mask_CNV &
#        ./muti_func_snp_compare.py -b 1000000 -i /data2/rawdata2/readDepth/YM05/chr${num}${i}.1M.norm --mask_cnv on -o /data/user/yangzz/mapping/09.filterVCF/GT_AD/MASK/lx99_mask/chr${num}${i}.mask_CNV &
#   #./muti_func_snp_compare.py -b 1000000 --high_threshold 8 --DP_low 3 --DP_high 20 --mask on -n ${sample1} --chromsome ${num}${i} -1 5 -i lx99_jm22/chr${num}${i}.gt_dp > MASK/lx99_mask/chr${num}${i}.mask_CNV_test &
#   #./muti_func_snp_compare.py -b 1000000 --high_threshold 6 --DP_low 3 --DP_high 20 --mask on -n ${sample2} --chromsome ${num}${i} -1 8 -i lx99_jm22/chr${num}${i}.gt_dp > MASK/jm22_mask/chr${num}${i}.mask_CNV_test &
#    done
#done

#for num in {1..7}; do
#    for i in {A,B,D};do
#        for a in {snp,indel}; do
#   path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
#            ./muti_func_snp_compare.py -b 1000000 -i lx99_jm22/chr${num}${i}_${a}_lxjm.gt --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 -2 8 --snp_graphy on --chromosome MASK/lx99_jm22_mask/chr${num}${i}.allCNV_use -o tmp/chr${num}${i}.lxjm_homo_${a}_density &
        #./muti_func_snp_compare.py -b 1000000 -i chr${num}${i}_snp_akcp.gt --DP_low 3 --DP_high 20 --GQ_sample 8 -1 5 -2 11 --snp_graphy on --chromosome CNV/chr${num}${i}.all_CNV -o akcp/chr${num}${i}.akcp01unmatch_homo_snp_density &
        #./muti_func_snp_compare.py -b 1000000 -i chr${num}${i}_snp_akcp.gt --DP_low 3 --DP_high 20 --GQ_sample 8 -1 11 -2 8 --snp_graphy on --chromosome CNV/chr${num}${i}.all_CNV -o akcp/chr${num}${i}.cp01cp03unmatch_homo_snp_density &
#        ./muti_func_snp_compare.py --density_graph_raw_vcf_two on -i chr${num}${i}_snp_akcp.gt -1 5 -2 8 --chromosome CNV/chr${num}${i}.all_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 --LEVEL "1.5,3.0" -o akcp/chr${num}${i}.akcp03unmatch_homo_snp_level &
#        ./muti_func_snp_compare.py --density_graph_raw_vcf_two on -i chr${num}${i}_snp_akcp.gt -1 5 -2 11 --chromosome CNV/chr${num}${i}.all_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 --LEVEL "1.5,3.0" -o akcp/chr${num}${i}.akcp01unmatch_homo_snp_level &
#        ./muti_func_snp_compare.py --density_graph_raw_vcf_two on -i chr${num}${i}_snp_akcp.gt -1 11 -2 8 --chromosome CNV/chr${num}${i}.all_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 --LEVEL "1.5,3.0" -o akcp/chr${num}${i}.cp01cp03unmatch_homo_snp_level &
#        done
#    done       
#done

#for num in {1..7}; do             # product snp density data
#    sample1=/data2/rawdata2/readDepth/YM05
#    sample2=/data2/rawdata2/readDepth/NZ10
#    path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
#    ./muti_func_snp_compare.py -b 1000000 --density on --DP_low 3 --DP_high 20 -n ${sample1} --chromsome ${num}A -1 5 -i lx99_jm22/chr${num}A.gt_dp > lx99_jm22/chr${num}A.lx99_density_ratio &

#    ./muti_func_snp_compare.py -b 1000000 --density on --DP_low 3 --DP_high 20 -n ${sample2} --chromsome ${num}A -1 8 -i lx99_jm22/chr${num}A.gt_dp > lx99_jm22/chr${num}A.jm22_density_ratio &
#done

#for num in {1..7}; do                                                 
#  path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
#  awk '{printf $1"\t"$2"\t"$2}{for (i=3;i<=NF;i++) {printf "\t"$i}printf "\n"}' lx99_jm22/chr${num}A.gt_dp > chr${num}A_lx99_jm22.bed & 
#done

#for num in {1..7}; do                                                 
#  path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
#  bedtools intersect -b chr${num}A_lx99_jm22.bed -a iwgsc_gene.bed -wb -sorted|awk '{printf $1"\t"}{for (i=7;i<=NF-1;i++){printf $i"\t"}{printf $NF}{printf "\n"}}'>chr${num}A_lx99_jm22_gene.gt & 
#done

#for num in {1..7}; do          # compare snp density detween two sample
#  ./muti_func_snp_compare.py -b 1000000 -i lx99_jm22/chr${num}A.gt_dp --DP_low "3,3" --DP_high "20,20" --LEVEL "0.4,0.8,1.0" -1 5 -2 8 --density_graphy on --chromsome ${num}A -o lx99_jm22/chr${num}A.lx99_density_graphy_gene & 

#  ./muti_func_snp_compare.py -b 1000000 -i lx99_jm22/chr${num}A.gt_dp --DP_low 3 --DP_high 20  -1 5 -2 8 --gene_snp_filter on --chromsome ${num}A -o chr${num}A.gene_pos &

#  ./muti_func_snp_compare.py -b 1000000 -i lx99_jm22/chr${num}A.gt_dp --DP_low "3,3" --DP_high "20,20" --centrosome "215,239,173,337,346,268,300,346,241,300,317,184,254,206,189,286,327,211,357,287,340" -1 5 -2 8 --long_short_count on --chromsome ${num}A >> test.test 
#  ./density_graph.py -b 1000000 -i chr${num}A_lx99_jm22.gt -1 5 --la 21185 --lb 11219 --lc 6065 -n ${num}A -o lx99_jm22/chr${num}A.jm22_density_graphy &
#done



#if [ ! -d "jm22_snp_indel_analysis/snp" ]; then
#     mkdir jm22_snp_indel_analysis/snp
#fi
#if [ ! -d "jm22_snp_indel_analysis/indel" ]; then
#     mkdir jm22_snp_indel_analysis/indel
#fi
#
#if [ ! -d "lx99_snp_indel_analysis/snp" ]; then
#     mkdir lx99_snp_indel_analysis/snp
#fi
#if [ ! -d "lx99_snp_indel_analysis/indel" ]; then
#     mkdir lx99_snp_indel_analysis/indel
#fi
#
#for num in {1..7}; do  # unmatch snp / delition snp number                                               
#    for i in {A,B,D};do
        #./muti_func_snp_compare.py --cnv_region on -i MASK/jm22_mask/chr${num}${i}.mask_CNV -o MASK/region/chr${num}${i}_jm22 &
        #./muti_func_snp_compare.py --cnv_region on -i MASK/lx99_mask/chr${num}${i}.mask_CNV -o MASK/region/chr${num}${i}_lx99 &
        #./muti_func_snp_compare.py --cnv_region on -i MASK/chr${num}${i}.bothCNV -o MASK/region/chr${num}${i}_both &
        #./muti_func_snp_compare.py --long_short_count on -i chr${num}${i}_snp_lxjm.gt --chromosome MASK/chr${num}${i}_all_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -o SNP/region/chr${num}${i}_lxjm_hohe &
        #./muti_func_snp_compare.py --long_short_count on -i chr${num}${i}_indel_lxjm.gt --chromosome MASK/chr${num}${i}_all_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -o Indel/region/chr${num}${i}_lxjm_hohe &
        #./muti_func_snp_compare.py --two_diff on -i chr${num}${i}_snp_lxjm.gt --chromosome MASK/chr${num}${i}_all_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -o SNP/region/chr${num}${i}_lxjm_diff &
        #./muti_func_snp_compare.py --two_diff on -i chr${num}${i}_indel_lxjm.gt --chromosome MASK/chr${num}${i}_all_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -o Indel/region/chr${num}${i}_lxjm_diff &
#        bcftools view -Ov /data/user/yangzz/mapping/S15/ReseqPipeNew/08.mergeGVCF/lx99_jm22/chr${num}${i}.ann.bcf.gz| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n'|./muti_func_snp_compare.py --per_snp_filter on  --chromosome MASK/lx99_jm22_mask/chr${num}${i}.allCNV_use --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -o lx99_jm22/region/annv1p1/chr${num}${i}.all_ann &
        #./muti_func_snp_compare.py --density_graph_raw_vcf_two on -i lx99_jm22/chr${num}${i}_snp_lxjm.gt --chromosome MASK/lx99_jm22_mask/chr${num}${i}.allCNV_use --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 --LEVEL "1.5,2.85" -o lx99_jm22/region/chr${num}${i}_snp_level_new &
        #cat jm22_snp_indel_analysis/snp/chr${num}${i}_snp_altLEVEL lx99_snp_indel_analysis/snp/chr${num}${i}_snp_altLEVEL|sort -n -k2|sed  '1,2d'|./filtlevel.py --filt  on --chromosome MASK/chr${num}${i}_all_CNV > lx99_jm22/region/chr${num}${i}_solo_level &
        #awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a)print a[$2]"\t"$4}' lx99_jm22/region/chr${num}${i}_solo_level lx99_jm22/region/chr${num}${i}_snp_level > lx99_jm22/region/chr${num}${i}_combine_solo_plus 
        #./filtlevel.py --compare on -i lx99_jm22/region/chr${num}${i}_combine_solo_plus -o lx99_jm22/region/chr${num}${i}_error
        #./muti_func_snp_compare.py --transition_transversion_count on -i chr${num}${i}_snp_lxjm.gt --chromosome lx99_jm22/region/chr${num}${i}_snp_level --LEVEL "1.5,2.85" --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 > lx99_jm22/region/chr${num}${i}_diff_transi_ver &
#        ./muti_func_snp_compare.py --indel_count on -i chr${num}${i}_indel_lxjm.gt --chromosome lx99_jm22/region/chr${num}${i}_snp_level --LEVEL "1.5,2.85" --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 | sort -nk2,2 > lx99_jm22/region/chr${num}${i}_jm22_indel_len_count &

#    done
#done
#wait
#cat lx99_jm22/region/chr*_error > lx99_jm22/region/combine_error
#cat lx99_jm22/region/chr*_snp_level > lx99_jm22/region/combine_snp_level
#cat SNP/region/chr*_lxjm_diff > SNP/region/combine_lxjm_diff_snp
#cat Indel/region/chr*_lxjm_diff > Indel/region/combine_lxjm_diff_indel

#cat  MASK/region/chr*_jm22 >  MASK/region/combine_jm22.region
#cat  MASK/region/chr*_lx99 >  MASK/region/combine_lx99.region
#cat  MASK/region/chr*_both >  MASK/region/combine_both.region

#for num in {1..7};do
#    for i in {A,B,D};do
        #for item in {snp,indel};do
#            echo "chr${num}${i}_${item}.gt"
#            ./muti_func_snp_compare.py -b 1000000 --density_graph_raw_vcf on \
#                                        --DP_low 3 --DP_high 20 --GQ_sample 8 \
#                                        --LEVEL 3 --sample /data2/rawdata2/readDepth/YM05 -1 5 -i lx99_snp_indel_analysis/chr${num}${i}_${item}.gt -o lx99_snp_indel_analysis/${item}/chr${num}${i}_${item}_altLEVEL &
#            echo "chr${num}${i}_${item}.gt"
#            ./muti_func_snp_compare.py -b 1000000 --density_graph_raw_vcf on \
#                                        --DP_low 3 --DP_high 20 --GQ_sample 8 \
#                                        --LEVEL 3 --sample /data2/rawdata2/readDepth/NZ10 -1 5 -i jm22_snp_indel_analysis/chr${num}${i}_${item}.gt -o jm22_snp_indel_analysis/${item}/chr${num}${i}_${item}_altLEVEL &


#            ./muti_func_snp_compare.py -b 1000000 --filter_and_count_snp on \
#                                       --DP_low 3 --DP_high 20 --GQ_sample 8 \
#                                       --sample /data2/rawdata2/readDepth/YM05 -1 5 -i lx99_snp_indel_analysis/chr${num}${i}_${item}.gt -o lx99_snp_indel_analysis/${item}/chr${num}${i}_${item}_altHOMO &
#            ./muti_func_snp_compare.py -b 1000000 --filter_and_count_snp on \
#                                       --DP_low 3 --DP_high 20 --GQ_sample 8 \
#                                       --sample /data2/rawdata2/readDepth/NZ10 -1 5 -i jm22_snp_indel_analysis/chr${num}${i}_${item}.gt -o jm22_snp_indel_analysis/${item}/chr${num}${i}_${item}_altHOMO &
#        item=indel
#        ./muti_func_snp_compare.py -b 1000000 --filter_and_count_indel on \
#                                   --DP_low 3 --DP_high 20 --GQ_sample 8 \
#                                   --sample /data2/rawdata2/readDepth/YM05 -1 5 -i lx99_snp_indel_analysis/chr${num}${i}_${item}.gt -o lx99_snp_indel_analysis/${item}/chr${num}${i}_${item}_altHOMO &
#        ./muti_func_snp_compare.py -b 1000000 --filter_and_count_indel on \
#                                   --DP_low 3 --DP_high 20 --GQ_sample 8 \
#                                   --sample /data2/rawdata2/readDepth/NZ10 -1 5 -i jm22_snp_indel_analysis/chr${num}${i}_${item}.gt -o jm22_snp_indel_analysis/${item}/chr${num}${i}_${item}_altHOMO &

        #done
        #wait_all
#    done
#    wait_all
#done
#wait
#cat lx99_snp_indel_analysis/snp/chr*_snp_altHOMO > lx99_snp_indel_analysis/snp/combine_snp_altHOMO
#cat lx99_snp_indel_analysis/indel/chr*_indel_altHOMO > lx99_snp_indel_analysis/indel/combine_indel_altHOMO
#cat jm22_snp_indel_analysis/snp/chr*_snp_altHOMO > jm22_snp_indel_analysis/snp/combine_snp_altHOMO
#cat jm22_snp_indel_analysis/indel/chr*_indel_altHOMO > jm22_snp_indel_analysis/indel/combine_indel_altHOMO

# ratio >= log10(-1.25)
#for num in {1..7};do
#    for i in {A,B,D};do
#        for item in {homo_DP,homo_GQ};do
#            ./muti_func_snp_filter.py -i chr${num}${i}_snp_akcp.gt --single_DP_GQ on -s 5 --chromsome CNV/Ak58Simu/chr${num}${i}.mask_CNV |grep ${item}|sort -nk1,1 > akcp/Ak58/chr${num}${i}_snp_${item} &
#
#            ./muti_func_snp_filter.py -i chr${num}${i}_snp_akcp.gt --single_DP_GQ on -s 8 --chromsome CNV/CP03/chr${num}${i}.mask_CNV |grep ${item}|sort -nk1,1> akcp/CP03/chr${num}${i}_snp_${item} &
#            ./muti_func_snp_filter.py -i chr${num}${i}_snp_akcp.gt --single_DP_GQ on -s 11 --chromsome CNV/CP01/chr${num}${i}.mask_CNV |grep ${item}|sort -nk1,1> akcp/CP01/chr${num}${i}_snp_${item} &
#
#            ./muti_func_snp_filter.py -i jm22_snp_indel_analysis/chr${num}${i}_indel.gt \
#                                      --single_DP_GQ on -s 5 --chromsome MASK/chr${num}A_combine_total_CNV \
#                                      |grep ${item} > jm22_snp_indel_analysis/indel/chr${num}${i}_indel_${item} &
#
#            ./muti_func_snp_filter.py -i lx99_snp_indel_analysis/chr${num}${i}_indel.gt \
#                                      --single_DP_GQ on -s 5 --chromsome MASK/chr${num}A_combine_total_CNV \
#                                      |grep ${item}> lx99_snp_indel_analysis/indel/chr${num}${i}_indel_${item} &
#            wait_all
#        done
#    done

#    ./muti_func_snp_filter.py -i lx99_jm22/chr${num}A_lx99_jm22_indel.gt> lx99_jm22/indel/chr${num}A_indel_filter.gt &
#    ./muti_func_snp_compare.py -b 1000000 -i lx99_jm22/chr${num}A.gt_dp --DP_low "3,3" --DP_high "20,20" --LEVEL "0.4,0.8,1.0" -1 5 -2 8 --density_graphy on --chromsome ${num}A -o lx99_jm22/chr${num}A.lx99_density_graphy_gene
#  ./lx99_jm22_snp_graphy.py -b 1000000 -i chr${num}A_lx99_jm22.gt -1 5 -2 7  -n ${num}A -o lx99_jm22/chr${num}A.density_filter &
#done
