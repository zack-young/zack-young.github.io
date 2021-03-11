CHRlst=chr1A.1,chr1A.2,chr1B.1,chr1B.2,chr1D.1,chr1D.2,chr2A.1,chr2A.2,chr2B.1,chr2B.2,chr2D.1,chr2D.2,chr3A.1,chr3A.2,chr3B.1,chr3B.2,chr3D.1,chr3D.2,chr4A.1,chr4A.2,chr4B.1,chr4B.2,chr4D.1,chr4D.2,chr5A.1,chr5A.2,chr5B.1,chr5B.2,chr5D.1,chr5D.2,chr6A.1,chr6A.2,chr6B.1,chr6B.2,chr6D.1,chr6D.2,chr7A.1,chr7A.2,chr7B.1,chr7B.2,chr7D.1,chr7D.2,chrUn.1,chrUn.2
#for CHR in `echo ${CHRlst} | tr "," " "`; do
#bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${CHR}.1k.bed -b ${CHR}.30_dedup.bam -counts -sorted > ${CHR}.1k.30_DP
#bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${CHR}.1k.bed -b ${CHR}.60_dedup.bam -counts -sorted > ${CHR}.1k.60_DP
#bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${CHR}.1k.bed -b ${CHR}.90_dedup.bam -counts -sorted > ${CHR}.1k.90_DP
#done

cat chr*.1M.30_tmp > combine_30_1M_DP
ave=`/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/count_CNV.py --normalize on -i combine_30_1M_DP`
for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${CHR}.1.1M.30_tmp > ${CHR}.1M.30_norm
  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${CHR}.2.1M.30_tmp >> ${CHR}.1M.30_norm
done
