#for num in {1..7}; do
#    sample1=/data2/rawdata2/readDepth/YM05
#    sample2=/data2/rawdata2/readDepth/NZ10
    #paste -d "\t" ${sample1}/chr${num}A.1M.norm ${sample2}/chr${num}A.1M.norm |sed '1d'|awk '{print $0"\t""chr""'$num'""A"}' > chr${num}A_combine.1M.norm
    #paste -d "\t" ${sample1}/chr${num}B.1M.norm ${sample2}/chr${num}B.1M.norm |sed '1d'|awk '{print $0"\t""chr""'$num'""B"}' > chr${num}B_combine.1M.norm
    #paste -d "\t" ${sample1}/chr${num}D.1M.norm ${sample2}/chr${num}D.1M.norm |sed '1d'|awk '{print $0"\t""chr""'$num'""D"}' > chr${num}D_combine.1M.norm
    #paste -d "\t" chr${num}A_combine.1M.norm ../lx99_jm22/chr${num}A.diff >> combine_gene_pos
    #paste -d "\t" chr${num}B_combine.1M.norm ../lx99_jm22/chr${num}B.diff >> combine_gene_pos
    #paste -d "\t" chr${num}D_combine.1M.norm ../lx99_jm22/chr${num}D.diff >> combine_gene_pos
#    grep "chr${num}A" test
#done

GFF=/data/annotation/wheat/CS_IWGSC/v1/commom/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3
#data=`cut -f9 ${GFF}|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u`
#echo `cut -f9 ${GFF}|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|wc -l`
for num in {1..7}; do
#    bedtools intersect -b gene_snp_indel.bed -a ../chr${num}A_lx99_jm22.bed -sorted > chr${num}A_snp.bed &
#    bedtools intersect -b gene_snp_indel.bed -a ../chr${num}B_lx99_jm22.bed -sorted > chr${num}B_snp.bed &
#    bedtools intersect -b gene_snp_indel.bed -a ../chr${num}D_lx99_jm22.bed -sorted > chr${num}D_snp.bed &
#    bedtools intersect -b gene_snp_indel.bed -a ../../chr${num}A.filter_vcf -sorted > chr${num}A_snp_alt_new.vcf &
#    bedtools intersect -b gene_snp_indel.bed -a ../../chr${num}B.filter_vcf -sorted > chr${num}B_snp_alt_new.vcf &
#    bedtools intersect -b gene_snp_indel.bed -a ../../chr${num}D.filter_vcf -sorted > chr${num}D_snp_alt_new.vcf &
#    cat vcf_header chr${num}A_snp_alt_new.vcf > chr${num}A_annotate.vcf
#    cat vcf_header chr${num}B_snp_alt_new.vcf > chr${num}B_annotate.vcf
#    cat vcf_header chr${num}D_snp_alt_new.vcf > chr${num}D_annotate.vcf

#    ../muti_func_snp_compare.py -i chr${num}A_snp.bed --DP_low "3,3" --DP_high "20,20"  -1 6 -2 9 --per_snp_filter on |grep -v ^$ > chr${num}A_snp_alt.bed &
#    ../muti_func_snp_compare.py -i chr${num}B_snp.bed --DP_low "3,3" --DP_high "20,20"  -1 6 -2 9 --per_snp_filter on |grep -v ^$ > chr${num}B_snp_alt.bed &
#    ../muti_func_snp_compare.py -i chr${num}D_snp.bed --DP_low "3,3" --DP_high "20,20"  -1 6 -2 9 --per_snp_filter on |grep -v ^$ > chr${num}D_snp_alt.bed &
#    bedtools intersect -b chr${num}A_snp_alt.bed -a ../../chr${num}A.bcf.gz -sorted > chr${num}A_snp_alt.vcf & 
#    bedtools intersect -b chr${num}B_snp_alt.bed -a ../../chr${num}B.bcf.gz -sorted > chr${num}B_snp_alt.vcf & 
#    bedtools intersect -b chr${num}D_snp_alt.bed -a ../../chr${num}D.bcf.gz -sorted > chr${num}D_snp_alt.vcf & 

#    java -Xmx8G -jar /data/user/yangzz/vcf/snpeff/snpEff/snpEff.jar -i vcf -o vcf tri1 chr${num}A_annotate.vcf > annotation/chr${num}A_snpeff.vcf
#    ./annotation/01CallTotalVariation.py -i annotation/chr${num}A_snpeff.vcf -o annotation/chr${num}A_snpeff_step1.txt
#    ./annotation/02callVatType.py annotation/chr${num}A_snpeff_step1.txt annotation/chr${num}A_snpeff_step2.txt
#
#     java -Xmx8G -jar /data/user/yangzz/vcf/snpeff/snpEff/snpEff.jar -i vcf -o vcf tri1 chr${num}D_annotate.vcf > annotation/chr${num}D_snpeff.vcf
#     ./annotation/01CallTotalVariation.py -i annotation/chr${num}D_snpeff.vcf -o annotation/chr${num}D_snpeff_step1.txt
#     ./annotation/02callVatType.py annotation/chr${num}D_snpeff_step1.txt annotation/chr${num}D_snpeff_step2.txt
#
#     java -Xmx8G -jar /data/user/yangzz/vcf/snpeff/snpEff/snpEff.jar -i vcf -o vcf tri1 chr${num}B_annotate.vcf > annotation/chr${num}B_snpeff.vcf
#     ./annotation/01CallTotalVariation.py -i annotation/chr${num}B_snpeff.vcf -o annotation/chr${num}B_snpeff_step1.txt
#     ./annotation/02callVatType.py annotation/chr${num}B_snpeff_step1.txt annotation/chr${num}B_snpeff_step2.txt

#    ./annotation/04asd.py annotation/chr${num}A_snpeff.vcf > annotation/chr${num}A_snpeff_step3.csv
#    echo "${num}A"
#    ./annotation/04asd.py annotation/chr${num}B_snpeff.vcf > annotation/chr${num}B_snpeff_step3.csv
#    echo "${num}B"
#    ./annotation/04asd.py annotation/chr${num}D_snpeff.vcf > annotation/chr${num}D_snpeff_step3.csv
#    echo "${num}D"
#    A1=`grep "${num}A" GO_use_level2|wc -l`
#    B1=`grep "${num}B" GO_use_level2|wc -l`
#    D1=`grep "${num}D" GO_use_level2|wc -l`
#
#    A2=`grep "${num}A" GO_use_level3|wc -l`
#    B2=`grep "${num}B" GO_use_level3|wc -l`
#    D2=`grep "${num}D" GO_use_level3|wc -l`
#    
#    echo "${num}Alevel2 ${A1} ${num}Alevel3 ${A2}" >> total_gene
#    echo "${num}Blevel2 ${B1} ${num}Blevel3 ${B2}" >> total_gene
#    echo "${num}Dlevel2 ${D1} ${num}Dlevel3 ${D2}" >> total_gene

done

