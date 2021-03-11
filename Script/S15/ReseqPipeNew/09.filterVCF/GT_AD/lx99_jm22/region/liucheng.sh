for i in {HC,LC}; do
#    for a in {jm22,lx99}; do
    for a in {mid,high}; do
        bedtools intersect -a /data/user/yangzz/mapping/09.filterVCF/GT_AD/basepart/IWGSC_v1.1_${i}_20170706_sorted.gff3 \
                           -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/combine_snp_${a} \
                           -sorted|grep 'gene'|sort -k9,9 -u |awk -F '[\t=;]' '{print$1"\t"$4"\t"$5"\t"$10}'> combine_snp_${a}_${i}.bed
        bedtools intersect -a /data/user/yangzz/mapping/09.filterVCF/GT_AD/basepart/IWGSC_v1.1_${i}_20170706_sorted.gff3 \
                           -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/combine_snp_${a} \
                           -sorted|cut -f9|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u > /data/user/yangzz/mapping/09.filterVCF/GT_AD/lx99_jm22/region/combine_snp_${a}_${i}gene
#        awk 'NR==FNR{a[$4]=$1"\t"$2"\t"$3;next}NR>FNR{if($1 in a) print a[$1]"\t"$0}' combine_snp_${a}_${i}.bed ~/IWGSC-ath-rice-v1p1.xls > combine_snp_${a}_${i}gene.ann

##    numB=`wc -l ${B}`
##    numD=`wc -l ${D}`
#    echo "${num}A ${A}" >> filter_gene
#    echo "${num}B ${B}" >> filter_gene
#    echo "${num}D ${D}" >> filter_gene
    done
done

#join -1 1 -2 1 -a1 gene_name IWGSC-ath-rice.ann > gene_name_annotate # find same part between files
#sed -i 's/ /\t/1' gene_name_annotate # substitute first blank into tab
#sort -k 1,1 -k 2n,2  GLM.bed > GLM.bednew  multiline sort
#bedtools intersect -a /data/user/yangzz/mapping/09.filterVCF/GT_AD/basepart/IWGSC_v1.1_HC_20170706_sorted.gff3 -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/MASK/lx99_jm22_mask/combine_jm22have -sorted|grep 'gene'|sort -k9,9 -u > jm22_CNV_HCgene.gff3
#awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a " " arr[a]}' combine_transi_ver_high
