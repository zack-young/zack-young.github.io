#!/bin/bash
set -euxo pipefail

while getopts 'l:' opt
do
    case "$opt" in
    l) meta_data=$OPTARG ;;
    esac
done
#time_date="2020-04-17162314" #"`date +%Y-%m-%d%H%M%S`" 
if [[ -z "${meta_data}" ]]; then
    echo "No first sample provided"; exit 1
fi
path_local="/data/user/yangzz/mapping/09.filterVCF/GT_AD"
arr=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
if true;then
for SAMPLE1 in ${arr[@]};do
    (if [[ ! -d "${SAMPLE1}" ]]; then
         mkdir ${SAMPLE1}
    fi
         for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do 
         #bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${chr}.1.1k.bed -b /data2/public/ReseqData/workflowFile/${SAMPLE1}/05.mergeAsplit/${chr}.*.bam -counts -sorted > ${SAMPLE1}/${chr}.1.1k.DP
         #bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${chr}.2.1k.bed -b /data2/public/ReseqData/workflowFile/${SAMPLE1}/05.mergeAsplit/${chr}.*.bam -counts -sorted > ${SAMPLE1}/${chr}.2.1k.DP

         gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' /data2/rawdata2/readDepth_ori_1k/${SAMPLE1}/${chr}.1.1k.DP > ${SAMPLE1}/${chr}.1.1M.tmp 
         gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' /data2/rawdata2/readDepth_ori_1k/${SAMPLE1}/${chr}.2.1k.DP > ${SAMPLE1}/${chr}.2.1M.tmp 
         done 
         cat ${SAMPLE1}/chr*.1M.tmp > ${SAMPLE1}/combine_1M_DP 
         
         sh ${path_local}/field_cultivar/script_1M.sh ${SAMPLE1}
    ) &
wait_all
done 
fi

for SAMPLE1 in ${arr[@]};do
      for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
     #1
      ${path_local}/muti_func_snp_compare.py -b 1000000 -i ${SAMPLE1}/${chr}.1M.norm --mask_cnv on -o ${SAMPLE1}/${chr}.mask_CNV &
     #2
      done
      wait 
      for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
      for item in {deletion,duplication};do
         grep "${item}" ${SAMPLE1}/${chr}.mask_CNV > ${SAMPLE1}/${chr}.mask_CNV_${item} &
      done
      done
wait_all
done
#sh bisample_compare_sequence_hugh_sample.sh -t ${time_date} -m ${meta_data} &
