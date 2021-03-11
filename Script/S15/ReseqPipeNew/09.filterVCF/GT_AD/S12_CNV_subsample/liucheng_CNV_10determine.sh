CHRlst=chr1A.1,chr1A.2,chr1B.1,chr1B.2,chr1D.1,chr1D.2,chr2A.1,chr2A.2,chr2B.1,chr2B.2,chr2D.1,chr2D.2,chr3A.1,chr3A.2,chr3B.1,chr3B.2,chr3D.1,chr3D.2,chr4A.1,chr4A.2,chr4B.1,chr4B.2,chr4D.1,chr4D.2,chr5A.1,chr5A.2,chr5B.1,chr5B.2,chr5D.1,chr5D.2,chr6A.1,chr6A.2,chr6B.1,chr6B.2,chr6D.1,chr6D.2,chr7A.1,chr7A.2,chr7B.1,chr7B.2,chr7D.1,chr7D.2,chrUn.1,chrUn.2
#for CHR in `echo ${CHRlst} | tr "," " "`; do
#bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${CHR}.1k.bed -b ${CHR}.10_dedup.bam -counts -sorted > ${CHR}.1k.10_DP
##bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${CHR}.1k.bed -b ${CHR}.60_dedup.bam -counts -sorted > ${CHR}.1k.60_DP
##bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${CHR}.1k.bed -b ${CHR}.90_dedup.bam -counts -sorted > ${CHR}.1k.90_DP
#done
#
#for chr in `echo ${CHRlst} | tr "," " "`;do
#gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' ${chr}.1k.10_DP > ${chr}.1M.10_tmp
#done
#
#cat chr*.1M.10_tmp > combine_10_1M_DP
#ave=`/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/count_CNV.py --normalize on -i combine_10_1M_DP`
#for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
#  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${CHR}.1.1M.10_tmp > ${CHR}.1M.10_norm
#  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${CHR}.2.1M.10_tmp >> ${CHR}.1M.10_norm
#  sed -i '1i 'S12'' ${CHR}.1M.10_norm
#done

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
    #1
#    sed -i '1i 'S12'' ${CHR}.1M.30_norm
#    sed -i '1i 'S12'' ${CHR}.1M.60_norm
#    sed -i '1i 'S12'' ${CHR}.1M.90_norm
    /data/user/yangzz/mapping/09.filterVCF/GT_AD/muti_func_snp_compare.py -b 1000000 -i ${CHR}.1M.10_norm --mask_cnv on -o ${CHR}.10_mask_CNV
   #2
    for item in {deletion,duplication};do
        grep "${item}" ${CHR}.10_mask_CNV > ${CHR}.10_mask_CNV_${item}
    done
done

