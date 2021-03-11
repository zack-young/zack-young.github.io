#!/bin/bash
set -euxo pipefail
while getopts 'c:v:a:d:m:t:M:' opt
do
    case "$opt" in
    c) SAMPLE1=$OPTARG ;;
    v) SAMPLE1_name=$OPTARG ;;
    a) SAMPLE1_vcf=$OPTARG ;;
    d) total=$OPTARG ;;
    m) meta_data=$OPTARG ;;
    t) time_date=$OPTARG ;;
    esac
done
if [[ -z "${SAMPLE1}" ]]; then
    echo "No first sample provided"; exit 1
fi
if [[ -z "${total}" ]]; then
    echo "No total sample bcf name provided."; exit 1
fi
if [[ -z "${meta_data}" ]]; then
    echo "No sample data provided."; exit 1
fi
if [[ -z "${time_date}" ]]; then
    echo "No time date provided."; exit 1
fi

SAMPLE_list=${SAMPLE1}
arr=(`echo ${SAMPLE_list} | tr ',' ' '`)

declare -A arr_2
while read line;do
arr_2[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $1}'`
done < ${meta_data}

declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $4}'`
done < ${meta_data}
#for num in {1..7};do
#    for i in {A,B,D};do
#        bcftools view  -v snps -Ov ~/mapping/08.mergeGVCF/${total}/chr${num}${i}.bcf.gz| bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > chr${num}${i}.snp_    matrix &
#    wait_all
#done
#done

#for sample1 in ${arr[@]};do


#for sample2 in ${arr_2[@]};do
#    unset arr_2[${sample1}]
#    if [ "${sample1}" != "${sample2}" ]; then


#while read line;do {
#    SAMPLE1=`echo $line | awk '{print $1}'`
#    num=0
#    for SAMPLE2 in ${arr[@]};do
#         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
#           if [[ ! -d "tmp/${SAMPLE1}_${SAMPLE2}_${time_date}" ]]; then
#                mkdir tmp/${SAMPLE1}_${SAMPLE2}_${time_date}
#           fi 
#           for num in {1..7}; do
#             for i in {A,B,D}; do
#               for item in {deletion,duplication};do
#                 a=`du -h tmp/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item}`
#                 b=`echo $a |awk '{print $1}'`
#                 if [ "${b}" == 0 ];then
#                   echo 'file zero' 
#                   awk '{print $1"\t"$2"\t"$3"\t"SAMPLE2"'"$item"'""_CNV"}' SAMPLE2="$SAMPLE2" tmp/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item}  > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item}
#                 else
#                   echo 'file not zero' 
#                   awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE2"'"$item"'""_CNV"}' SAMPLE2="$SAMPLE2" tmp/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item} tmp/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item}
#                 fi
#        #
#                 a=`du -h tmp/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item}`
#                 b=`echo $a |awk '{print $1}'`
#                 echo ${b}
#                 if [ "${b}" == 0 ];then
#                    echo 'file zero'
#                    awk '{print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" tmp/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item}  > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
#                 else
#                    echo 'file not zero'
#                    awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" tmp/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item} tmp/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
#                 fi
#                 awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==1) print $1"\t"$2"\t"$3"\t""'"$item"'""_both_CNV"}' tmp/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item} tmp/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}both_CNV_${item}
#            #3  
#                 cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item} tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item} tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}both_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV_${item}
#               done
#               cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV_deletion tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV_duplication | sort -u -nk2,2 > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV
#             done
#           done
#         else
#         unset arr[${num}]
#         arr=(`echo ${arr[@]}`)
#         fi
#         num=$(($num + 1))
##
#    done
#} &
#wait_all
#done < ${meta_data}
#wait
#
#for f_num in `seq 1 10 | tr '\n' ' '|sed s/,$//g`;do
#for num in {1..7};do
#    for i in {A,B,D};do
#        #for item in {snp,indel};do
#        ./Huge_samples_snp_compare.py -i chr${num}${i}.snp_matrix_${f_num} --single_DP_GQ on --meta_data ${meta_data} --time_data ${time_data} --suffix _${f_num} --chromosome chr${num}${i}  &
#        wait_all
#    done
#done
#done
#wait

## level "1.5,2.85"
#for f_num in `seq 1 20 | tr '\n' ' '|sed s/,$//g`;do
#for num in {1..7}; do
#    for i in {A,B,D};do
#        ./Huge_samples_snp_compare.py -b 1000000 -i chr${num}${i}.snp_matrix_${f_num} --meta_data ${meta_data}  --DP_low 3 --DP_high 99 --GQ_sample 8 --time_data ${time_date} --sample ${SAMPLE_list} --snp_graphy on --chromosome chr${num}${i} --suffix _${f_num} &
#        wait_all
#    done
#done
#done


for f_num in `seq 15 20 | tr '\n' ' '|sed s/,$//g`;do
for num in {1..7}; do
    for i in {A,B,D};do
        ./Huge_samples_snp_compare.py --density_graph_raw_vcf_two on -i ${total}/chr${num}${i}.snp_matrix_${f_num} --meta_data ${meta_data} --total ${total} --sample ${SAMPLE_list} --time_data     ${time_date} --chromosome chr${num}${i} --DP_low 3 --DP_high 99 --GQ_sample 8 -b 1000000 --LEVEL "19" --suffix _${f_num} &
        wait_all
    done
done
done
wait
#cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr*.${SAMPLE1}_${SAMPLE2}_unmatchhomo_snp_density |grep -v 'diff' > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/combine.${SAMPLE1}_${SAMPLE2}_unm    atchhomo_snp_density


#for sample1 in ${arr[@]};do


#for sample2 in ${arr_2[@]};do
#    unset arr_2[${sample1}]
#    if [ "${sample1}" != "${sample2}" ]; then
#    for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do

    #cat tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density_* > tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density
    #./muti_func_snp_compare.py --two_diff_level on -i tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density --sample1 ${sample1} --sample2 ${sample2} --time_data ${time_date} --total ${total} &
#    a1=`du -h tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.${sample1}to${sample2}all_CNV`
#    b1=`echo ${a1} |awk '{print $1}'`
    #cat tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.only_homo_snp_level tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.${sample1}to${sample2}all_CNV > tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_level
    #if [ "${b}" != 0 ];then
    #awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" tmp/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item} tmp/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
    #fi
   #awk '{print $4}' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_density > tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_density
   #sed -i '1i '${sample}_${SAMPLE2}'' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_density
#    done
#    ./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${total}/${sample1}_${sample2}_${time_date}" -s ".homo_snp_level" > plotfile/plotfile_${sample1}_${sample2}
#    Rscript ~/R/graph_diff_sim.R --sample1 ${sample1} --sample1_name ${arr_3[${sample1}]} --sample2 ${sample2} --sample2_name ${arr_3[${sample2}]} -s ${time_date}
#    fi
#paste `awk '{if ($1!="'"${sample}"'")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/""'"${sample}"'""_"$1"_2019-12-15201714/combine_homodiff_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_${sample}_${SAMPLE2}_homodiff_snp_density
#done
#cat `awk '{if ($1!="PH09")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/PH09_"$1"_2019-12-15201714/combine_homo_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_PH09_homo_snp_density
#done




#for sample in ${arr[@]};do
#while read line;do
#    SAMPLE2=`echo $line | awk '{print $1}'`
#    if [ "${sample}" != "${SAMPLE2}" ]; then
#    if [[ -s "tmp/${sample}_${SAMPLE2}_${time_date}/chr1A.homo_snp_level_20" ]]; then
#    cat tmp/${sample}_${SAMPLE2}_${time_date}/chr*.homo_snp_level_* > tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_level
#    awk '{print $5}' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_level > tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_level
#    sed -i '1i '${sample}_${SAMPLE2}'' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_level
#    cp tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_level field_cultivar/tmp_snp_level/${sample}_${SAMPLE2}_homo_snp_level
#    cp tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_level field_cultivar/tmp_snp_level/${sample}_${SAMPLE2}_homodiff_snp_level
#    fi
#    fi
#done < ${meta_data}
#paste `awk '{if ($1!="'"${sample}"'")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/""'"${sample}"'""_"$1"_2019-12-15201714/combine_homodiff_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_${sample}_${SAMPLE2}_homodiff_snp_density
#done
#cat `awk '{if ($1!="PH09")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/PH09_"$1"_2019-12-15201714/combine_homo_snp_density"}' field_cultivar/    metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_PH09_homo_snp_density



#while read line;do
#a=`echo $line | awk '{print $1}'`
#b=`echo $line | awk '{print $4}'` 
#./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${SAMPLE1}_${a}_${time_date}" -s ".${SAMPLE1}_${b}unmatch_homo_snp_level" > plotfile/plotfile_${SAMPLE1}_${b}
#./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${SAMPLE1}_${a}_${time_date}" -s ".${SAMPLE1}to${a}all_CNV_deletion" > plotfile/plotfile_${SAMPLE1}_${a}_deletion
#./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${SAMPLE1}_${a}_${time_date}" -s ".${SAMPLE1}to${a}all_CNV_duplication" > plotfile/plotfile_${SAMPLE1}_${a}_duplication
#Rscript ~/R/graph_diff.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${a} --sample2_name ${b} -s ${time_date}
#Rscript ~/R/graph_CNV_deletion.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${a} --sample2_name ${b} -s ${time_date}
#Rscript ~/R/graph_CNV_duplication.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${a} --sample2_name ${b} -s ${time_date}
#done < field_cultivar/${meta_data}
#Rscript ~/R/graph_CNV.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${SAMPLE2} --sample2_name ${SAMPLE2_name} -s ${time_date}
#for num in {1..7}; do
#    for i in {A,B,D};do
#        ./muti_func_snp_compare.py --transition_transversion_count on -i ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_snp.gt --chromosome tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}.${SAMPLE1}_${SAMPLE2}unmatch_homo_snp_level  --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -1 ${first} -2 ${second} | sort -nk2,2 > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_diff_transi_ver &
#        ./muti_func_snp_compare.py --indel_count on -i ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_indel.gt --chromosome tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}.${SAMPLE1}_${SAMPLE2}unmatch_homo_snp_level --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -1 ${first} -2 ${second}  | sort -nk2,2 > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_indel_len_count &
#    done
#done
#wait
#cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr*_diff_transi_ver | awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/combine_total_diff_transi_ver
#cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr*_diff_transi_ver |grep -E 'high|mid'| awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/combine_diff_diff_transi_ver
#cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr*_indel_len_count | awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/combine_total_indel_len_count 
#cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr*_indel_len_count  |grep -E 'high|mid'| awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/combine_diff_indel_len_count 
#
#Rscript ~/R/transi_tranver.R -s ${SAMPLE1}_${SAMPLE2}_${time_date}
#Rscript ~/R/indel_length.R -s ${SAMPLE1}_${SAMPLE2}_${time_date}
