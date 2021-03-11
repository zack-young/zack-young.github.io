#!/bin/bash
for num in {1..7};do
#    cat ../MASK/chr${num}A_combine_total chr${num}A.lx99_density_graphy_gene|sort -nk2> chr${num}A_graphy
#    cat ../MASK/chr${num}B_combine_total chr${num}B.lx99_density_graphy_gene|sort -nk2> chr${num}B_graphy
#    cat ../MASK/chr${num}D_combine_total chr${num}D.lx99_density_graphy_gene|sort -nk2> chr${num}D_graphy 
#    grep 'level2' ../chr${num}A.lx99_density_graphy_gene > chr${num}A_level2_bin.bed
#    grep 'level2' ../chr${num}B.lx99_density_graphy_gene > chr${num}B_level2_bin.bed
#    grep 'level2' ../chr${num}D.lx99_density_graphy_gene > chr${num}D_level2_bin.bed
#    grep 'level3' ../chr${num}A.lx99_density_graphy_gene > chr${num}A_level3_bin.bed
#    grep 'level3' ../chr${num}B.lx99_density_graphy_gene > chr${num}B_level3_bin.bed
#    grep 'level3' ../chr${num}D.lx99_density_graphy_gene > chr${num}D_level3_bin.bed
    grep 'DP' ../chr${num}A_graphy > chr${num}A_DP_bin.bed
    grep 'DP' ../chr${num}B_graphy > chr${num}B_DP_bin.bed
    grep 'DP' ../chr${num}D_graphy > chr${num}D_DP_bin.bed
    grep 'unmatchCNV' ../chr${num}A_graphy > chr${num}A_CNV_bin.bed
    grep 'unmatchCNV' ../chr${num}B_graphy > chr${num}B_CNV_bin.bed
    grep 'unmatchCNV' ../chr${num}D_graphy > chr${num}D_CNV_bin.bed
    
done

#for num in {1..7};do
#    grep -F -f lx99_mask/chr${num}A.mask_CNV_test jm22_mask/chr${num}A.mask_CNV_test | sort -n -k3 | uniq > chr${num}A.mask_filtered_combine_test
#    grep -F -f lx99_mask/chr${num}B.mask_CNV_test jm22_mask/chr${num}B.mask_CNV_test | sort -n -k3 | uniq > chr${num}B.mask_filtered_combine_test
#    grep -F -f lx99_mask/chr${num}D.mask_CNV_test jm22_mask/chr${num}D.mask_CNV_test | sort -n -k3 | uniq > chr${num}D.mask_filtered_combine_test
#done

GFF=/data/annotation/wheat/CS_IWGSC/v1/commom/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3
#data=`cut -f9 ${GFF}|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u`
#echo `cut -f9 ${GFF}|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|wc -l`
for num in {1..7}; do
#    bedtools intersect -b chr${num}A_level2_bin.bed -a ${GFF} -sorted > gene2A 
#    bedtools intersect -b chr${num}B_level2_bin.bed -a ${GFF} -sorted > gene2B 
#    bedtools intersect -b chr${num}D_level2_bin.bed -a ${GFF} -sorted > gene2D 
#
#    bedtools intersect -b chr${num}A_level3_bin.bed -a ${GFF} -sorted > gene3A
#    bedtools intersect -b chr${num}B_level3_bin.bed -a ${GFF} -sorted > gene3B
#    bedtools intersect -b chr${num}D_level3_bin.bed -a ${GFF} -sorted > gene3D

    bedtools intersect -b chr${num}A_DP_bin.bed -a ${GFF} -sorted > gene2A
    bedtools intersect -b chr${num}B_DP_bin.bed -a ${GFF} -sorted > gene2B
    bedtools intersect -b chr${num}D_DP_bin.bed -a ${GFF} -sorted > gene2D
   
    bedtools intersect -b chr${num}A_CNV_bin.bed -a ${GFF} -sorted > gene3A
    bedtools intersect -b chr${num}B_CNV_bin.bed -a ${GFF} -sorted > gene3B
    bedtools intersect -b chr${num}D_CNV_bin.bed -a ${GFF} -sorted > gene3D


    A2=`cut -f9 gene2A|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|grep "${num}A"|wc -l` 
    B2=`cut -f9 gene2B|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|grep "${num}B"|wc -l`
    D2=`cut -f9 gene2D|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|grep "${num}D"|wc -l`

    A3=`cut -f9 gene3A|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|grep "${num}A"|wc -l`
    B3=`cut -f9 gene3B|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|grep "${num}B"|wc -l`
    D3=`cut -f9 gene3D|cut -d";" -f1|cut -d"=" -f2|sed 's/\..*//g'|sort -k1,1 -u|grep "${num}D"|wc -l`
#    D=`grep ${num}'D' ${data}`
#    numA=`wc -l ${A}`
#    numB=`wc -l ${B}`
#    numD=`wc -l ${D}`
    echo "${num}A ${A2}" >> DP_gene
    echo "${num}B ${B2}" >> DP_gene
    echo "${num}D ${D2}" >> DP_gene
    echo "${num}A ${A3}" >> CNV_gene
    echo "${num}B ${B3}" >> CNV_gene
    echo "${num}D ${D3}" >> CNV_gene
done

