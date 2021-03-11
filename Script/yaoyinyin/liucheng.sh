#!/usr/bin/env bash
#for num in {3,4};do
#num=3
#for i in {B,D}; do
#    bcftools view /data/user/yangzz/mapping/08.mergeGVCF/chr${num}${i}.ann.bcf.gz|bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t[%SAMPLE\t%GT\t]\n' > chr${num}${i}.ann.all &
#        wait_all
#    awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a == 1) print a[$2]"\t"$4}' ../jm22_snp_indel_analysis/snp/chr${num}${i}_snp_altLEVEL ../lx99_snp_indel_analysis/snp/chr${num}${i}_snp_altLEVEL|sed '1d' > reg    ion/chr${num}${i}_solo_level_1.1
    #awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a == 1) print a[$2]"\t"$4}' region/chr${num}${i}_solo_level_1.1 region/chr${num}${i}_snp_level  > region/chr${num}${i}_combine_solo_plus_1.1
#     sed -i "1iCHR     start   end     alt_ratio_level_jm22    alt_ratio_level_lx99" region/chr${num}${i}_solo_level_1.1 
#     sed -i "1iCHR     start   end     alt_ratio_level_jm22    alt_ratio_level_lx99    diff_homosnp_ratio_level" region/chr${num}${i}_combine_solo_plus_1.1

    #../filtlevel.py --compare on -i region/chr${num}${i}_combine_solo_plus_1.1 >  region/chr${num}${i}_count_solo_plus

#done


#while read line;do
#    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t[%GT\t%AD\t%GQ\t]\n' /data/user/yangzz/mapping/08.mergeGVCF/henan_shandong/chr1A.ann.bcf.gz| grep ${line} |grep -v 'inter' >> 1A_gene.vcf 
#done < 1A_gene.list
#while read line;do
#    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t[%GT\t%AD\t%GQ\t]\n' /data/user/yangzz/mapping/08.mergeGVCF/henan_shandong/chr1B.ann.bcf.gz| grep ${line} |grep -v 'inter' >> 1B_gene.vcf
#done < 1B_gene.list
#while read line;do
#    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t[%GT\t%AD\t%GQ\t]\n' /data/user/yangzz/mapping/08.mergeGVCF/henan_shandong/chr1D.ann.bcf.gz| grep ${line} |grep -v 'inter' >> 1D_gene.vcf
#done < 1D_gene.list
#while read line;do
#    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t[%GT\t%AD\t%GQ\t]\n' /data/user/yangzz/mapping/08.mergeGVCF/henan_shandong/chr6A.ann.bcf.gz| grep ${line} |grep -v 'inter' >> 6A_gene.vcf
#done < 6A_gene.list
#while read line;do
#    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t[%GT\t%AD\t%GQ\t]\n' /data/user/yangzz/mapping/08.mergeGVCF/henan_shandong/chr6B.ann.bcf.gz| grep ${line} |grep -v 'inter' >> 6B_gene.vcf
#done < 6B_gene.list
#while read line;do
#    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t[%GT\t%AD\t%GQ\t]\n' /data/user/yangzz/mapping/08.mergeGVCF/henan_shandong/chrUn.ann.bcf.gz| grep ${line} |grep -v 'inter' >> Un_gene.vcf
#done < Un_gene.list
sed -n '1p' total.gene.ann.final.vcf > HMW_gene_final
while read line;do
    cat total.gene.ann.final.vcf|grep ${line} >> HMW_gene_final
done < HMW.list
sed -n '1p' total.gene.ann.final.vcf > LMW_gene_final
while read line;do
    cat total.gene.ann.final.vcf|grep ${line} >> LMW_gene_final
done < LMW.list
sed -n '1p' total.gene.ann.final.vcf >  prolamin_gene_final
while read line;do
    cat total.gene.ann.final.vcf|grep ${line} >> prolamin_gene_final
done < prolamin.list
