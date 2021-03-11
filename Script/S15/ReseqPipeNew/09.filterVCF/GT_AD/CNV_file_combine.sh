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


for sample1 in ${arr[@]};do
for sample2 in ${arr_2[@]};do
    unset arr_2[${sample1}]
    if [ "${sample1}" != "${sample2}" ]; then

    (cat tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density_* > tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density
    ./muti_func_snp_compare.py --two_diff_level on -i tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density --sample1 ${sample1} --sample2 ${sample2} --time_data ${time_date} --total ${total}
    a1=`du -h tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.${sample1}to${sample2}all_CNV`
    b1=`echo ${a1} |awk '{print $1}'`
    cat tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.only_homo_snp_level tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.${sample1}to${sample2}all_CNV > tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_level) &
    wait_all
  #awk '{print $4}' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_density > tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_density
  #sed -i '1i '${sample}_${SAMPLE2}'' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_density
    wait
    ./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${total}/${sample1}_${sample2}_${time_date}" -s ".homo_snp_level" > plotfile/plotfile_${sample1}_${sample2}
    Rscript ~/R/graph_diff_sim.R --sample1 ${sample1} --sample1_name ${arr_3[${sample1}]} --sample2 ${sample2} --sample2_name ${arr_3[${sample2}]} -s ${time_date}
    fi
#paste `awk '{if ($1!="'"${sample}"'")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/""'"${sample}"'""_"$1"_2019-12-15201714/combine_homodiff_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_${sample}_${SAMPLE2}_homodiff_snp_density
done
#cat `awk '{if ($1!="PH09")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/PH09_"$1"_2019-12-15201714/combine_homo_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_PH09_homo_snp_density
done

