for i in NZ09 NZ28 both; do
    for a in HC LC; do
        grep "${i}_CNV" chr2D.NZ09_NZ28unmatch_homo_snp_level | bedtools intersect -a /data/user/yangzz/mapping/09.filterVCF/GT_AD/basepart/IWGSC_v1.1_${a}_20170706_sorted.gff3 -b - -sorted|cut -f9|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u > chr2D_${i}CNV_${a}_gene &
        grep "${i}_CNV" chr7D.NZ09_NZ28unmatch_homo_snp_level | bedtools intersect -a /data/user/yangzz/mapping/09.filterVCF/GT_AD/basepart/IWGSC_v1.1_${a}_20170706_sorted.gff3 -b - -sorted|cut -f9|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u > chr7D_${i}CNV_${a}_gene &
    done
done
