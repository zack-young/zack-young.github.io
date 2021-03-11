#!/bin/bash
set -x

while getopts 'c:v:a:s:m:b:d:f:n:l:' opt
do
    case "$opt" in
    c) SAMPLE=$OPTARG ;;
    d) total=$OPTARG ;;
    l) meta_data=$OPTARG ;;
    esac
done
time_date="`date +%Y-%m-%d%H%M%S`" 
if [[ -z "${SAMPLE}" ]]; then
    echo "No first sample provided"; exit 1
fi
if [[ -z "${total}" ]]; then
    echo "No total sample bcf name provided."; exit 1
fi
if [[ ! -d "tmp" ]]; then
     mkdir tmp
fi

SAMPLE_list=${SAMPLE}
arr=(`echo ${SAMPLE_list} | tr ',' ' '`)
for SAMPLE1 in ${arr[@]};do
    (if [[ ! -d "tmp/${total}/${SAMPLE1}_${time_date}" ]]; then
         mkdir tmp/${total}/${SAMPLE1}_${time_date}
    fi

    if [[ ! -d "field_cultivar/${SAMPLE1}" ]]&&[[ ! -d "CNV_heatmap/${SAMPLE1}_${time_date}" ]]; then
         mkdir CNV_heatmap/${SAMPLE1}_${time_date}
         path_SM="/data2/rawdata2/readDepth_ori_1k/${SAMPLE1}"
         # norm
         for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do 
         gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' $path_SM/${chr}.1.1k.DP > CNV_heatmap/${SAMPLE1}_${time_date}/${chr}.1.1M.tmp 
         gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' $path_SM/${chr}.2.1k.DP > CNV_heatmap/${SAMPLE1}_${time_date}/${chr}.2.1M.tmp 
         done 
         cat CNV_heatmap/${SAMPLE1}_${time_date}/chr*.1M.tmp > CNV_heatmap/${SAMPLE1}_${time_date}/combine_1M_DP 
         sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/script_1M.sh CNV_heatmap/${SAMPLE1}_${time_date} &
    fi 
    wait

    for num in {1..7}; do
        for i in {A,B,D}; do
        #1
        if [[ -d "field_cultivar/${SAMPLE1}" ]]&&[[ -s "field_cultivar/${SAMPLE1}/chr${num}${i}.1M.norm" ]]; then
            ./muti_func_snp_compare.py -b 1000000 -i field_cultivar/${SAMPLE1}/chr${num}${i}.1M.norm --mask_cnv on -o tmp/${total}/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV 
        else
            ./muti_func_snp_compare.py -b 1000000 -i CNV_heatmap/${SAMPLE1}_${time_date}/chr${num}${i}.1M.norm --mask_cnv on -o tmp/${total}/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV
        fi
        #2
        for item in {deletion,duplication};do
            grep "${item}" tmp/${total}/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV > tmp/${total}/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item}
        done
        done
    done) &
wait_all
done 

wait
#sh bisample_compare_sequence_hugh_sample.sh -t ${time_date} -m ${meta_data} &
