#!/bin/bash
set -euxo pipefail

while getopts 'l:' opt
do
    case "$opt" in
    l) meta_data=$OPTARG ;;
    esac
done
time_date="2020-04-17162314" #"`date +%Y-%m-%d%H%M%S`" 
if [[ -z "${meta_data}" ]]; then
    echo "No first sample provided"; exit 1
fi
path_local="/data/user/yangzz/mapping/09.filterVCF/GT_AD"
arr=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
if true;then
for SAMPLE1 in ${arr[@]};do
    (if [[ ! -d "${SAMPLE1}_${time_date}" ]]; then
         mkdir ${SAMPLE1}_${time_date}
    fi
         #for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do 
         #bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${chr}.1.1k.bed -b /data2/public/ReseqData/workflowFile/${SAMPLE1}/05.mergeAsplit/${chr}.*.bam -counts -sorted > ${SAMPLE1}_${time_date}/${chr}.1.1k.DP
         #bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${chr}.2.1k.bed -b /data2/public/ReseqData/workflowFile/${SAMPLE1}/05.mergeAsplit/${chr}.*.bam -counts -sorted > ${SAMPLE1}_${time_date}/${chr}.2.1k.DP
         #done 
         cat ${SAMPLE1}_${time_date}/chr*.1k.DP|cut -f4 > ${SAMPLE1}_${time_date}/combine_1k_DP 
         sh script_1k.sh ${SAMPLE1}_${time_date}
    ) &
wait_all
done 
fi

#sh bisample_compare_sequence_hugh_sample.sh -t ${time_date} -m ${meta_data} &
