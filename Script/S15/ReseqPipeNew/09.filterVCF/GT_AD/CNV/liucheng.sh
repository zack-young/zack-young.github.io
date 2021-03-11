for num in {1..7}; do
    for i in {A,B,D}; do 
#    sample1=/data2/rawdata2/readDepth/YM05
#    sample2=/data2/rawdata2/readDepth/NZ10
        #1
        #paste -d "\t" ../CNV_heatmap/CP01/chr${num}${i}.1M.norm ../CNV_heatmap/CP03/chr${num}${i}.1M.norm ../CNV_heatmap/Ak58Simu/chr${num}${i}.1M.norm|awk '{print $0"\t""chr""'$num'""'$i'"}' > chr${num}${i}_combine.1M.norm
        #2
        #../muti_func_snp_compare.py -b 1000000 -i ../CNV_heatmap/CP01/chr${num}${i}.1M.norm --mask_cnv on -o CP01/chr${num}${i}.mask_CNV 
        #../muti_func_snp_compare.py -b 1000000 -i ../CNV_heatmap/CP03/chr${num}${i}.1M.norm --mask_cnv on -o CP03/chr${num}${i}.mask_CNV 
        #../muti_func_snp_compare.py -b 1000000 -i ../CNV_heatmap/Ak58Simu/chr${num}${i}.1M.norm --mask_cnv on -o Ak58Simu/chr${num}${i}.mask_CNV 
        #3
        #awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $0"\t""lx987_CNV"}' 3097/chr${num}${i}.mask_CNV lx987/chr${num}${i}.mask_CNV > chr${num}${i}_lx987to3097_own_CNV
        #awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $0"\t""3097_CNV"}' lx987/chr${num}${i}.mask_CNV 3097/chr${num}${i}.mask_CNV > chr${num}${i}_3097tolx987_own_CNV
        #cat lx987/chr${num}${i}.mask_CNV 3097/chr${num}${i}.mask_CNV|awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $0"\t""5181_CNV"}' - 5181/chr${num}${i}.mask_CNV > chr${num}${i}_5181to30lx_own_CNV
        #4
        cat CP03/chr${num}${i}.mask_CNV CP01/chr${num}${i}.mask_CNV Ak58Simu/chr${num}${i}.mask_CNV | sort -u -nk2,2 > chr${num}${i}.all_CNV
        #paste -d "\t" ${sample1}/chr${num}${i}.1M.norm ${sample2}/chr${num}${i}.1M.norm |sed '1d'|awk '{print $0"\t""chr""'$num'""'$i'"}' > chr${num}${i}_combine.1M.norm
        #paste -d "\t" chr${num}${i}_combine.1M.norm ../lx99_jm22/chr${num}${i}.diff >> combine_gene_pos
   # grep "chr${num}A" test
    done
done

#GFF=/data/annotation/wheat/CS_IWGSC/v1/commom/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3

#data=`cut -f9 ${GFF}|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u`
#echo `cut -f9 ${GFF}|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|wc -l`


#for num in {1..7}; do
#    A=`cut -f9 ${GFF}|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|grep "${num}A"|wc -l` 
##    D=`grep ${num}'D' ${data}`
##    numA=`wc -l ${A}`
##    numB=`wc -l ${B}`
##    numD=`wc -l ${D}`
#    echo "${num}A ${A}" >> total_gene
#    echo "${num}B ${B}" >> total_gene
#    echo "${num}D ${D}" >> total_gene
#done
#
#gene=combine_gene

#bedtools intersect -a /data/user/yangzz/mapping/09.filterVCF/GT_AD/basepart/IWGSC_v1.1_HC_20170706_sorted.gff3 \
#                   -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/MASK/lx99_jm22_mask/combine_lx99have \
#                   -sorted|cut -f9|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u > gene/lx99_CNV_HCgene
#bedtools intersect -a /data/user/yangzz/mapping/09.filterVCF/GT_AD/basepart/IWGSC_v1.1_HC_20170706_sorted.gff3 \
#                   -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/MASK/lx99_jm22_mask/combine_jm22have \
#                   -sorted|cut -f9|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u > gene/jm22_CNV_HCgene
#
#bedtools intersect -a /data/user/yangzz/mapping/09.filterVCF/GT_AD/basepart/IWGSC_v1.1_LC_20170706_sorted.gff3 \
#                   -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/MASK/lx99_jm22_mask/combine_lx99have \
#                   -sorted|cut -f9|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u > gene/lx99_CNV_LCgene 
#bedtools intersect -a /data/user/yangzz/mapping/09.filterVCF/GT_AD/basepart/IWGSC_v1.1_LC_20170706_sorted.gff3 \
#                   -b /data/user/yangzz/mapping/09.filterVCF/GT_AD/MASK/lx99_jm22_mask/combine_jm22have \
#                   -sorted|cut -f9|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u > gene/jm22_CNV_LCgene






#for num in {1..7}; do
#    ./wanxin.py chr${num}A_combine.1M.norm chr${num}A_combine_unmatch.1M.norm chr${num}A_combine_own.1M.norm
#    cat chr${num}A_combine_own.1M.norm chr${num}A_combine_unmatch.1M.norm |sort -n -k1,1|awk '{print$5"\t"$1"\t"$2"\t"$6}'> chr${num}A_combine_total_CNV
#    A=`cut -f9 ${gene}|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|grep "${num}A"|wc -l`
#    B=`cut -f9 ${gene}|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|grep "${num}B"|wc -l`
#    D=`cut -f9 ${gene}|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|grep "${num}D"|wc -l`
##    D=`grep ${num}'D' ${data}`
##    numA=`wc -l ${A}`
##    numB=`wc -l ${B}`
##    numD=`wc -l ${D}`
#    echo "${num}A ${A}" >> filter_gene
#    echo "${num}B ${B}" >> filter_gene
#    echo "${num}D ${D}" >> filter_gene
#done

#join -1 1 -2 1 -a1 gene_name IWGSC-ath-rice.ann > gene_name_annotate # find same part between files
#sed -i 's/ /\t/1' gene_name_annotate # substitute first blank into tab
#sort -k 1,1 -k 2n,2  GLM.bed > GLM.bednew  multiline sort
