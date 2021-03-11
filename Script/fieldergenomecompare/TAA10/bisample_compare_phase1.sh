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
declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
done < ${meta_data}


if true;then
for SAMPLE1 in ${arr[@]};do
    (if [[ ! -d "${SAMPLE1}_${time_date}" ]]; then
         mkdir ${SAMPLE1}_${time_date}
     fi
         for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do 
         bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${chr}.1.1M.bed -b /data/user/yangzz/mapping/Reseq_data/${SAMPLE1}/05.mergeAsplit/${chr}.1.*.bam -counts -sorted > ${SAMPLE1}_${time_date}/${chr}.1.1M.DP
         bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${chr}.2.1M.bed -b /data/user/yangzz/mapping/Reseq_data/${SAMPLE1}/05.mergeAsplit/${chr}.2.*.bam -counts -sorted > ${SAMPLE1}_${time_date}/${chr}.2.1M.DP
         gawk '{print $4*1000/($3-$2)}' ${SAMPLE1}_${time_date}/${chr}.1.1M.DP > ${SAMPLE1}_${time_date}/${chr}.1.1M.tmp 
         gawk '{print $4*1000/($3-$2)}' ${SAMPLE1}_${time_date}/${chr}.2.1M.DP > ${SAMPLE1}_${time_date}/${chr}.2.1M.tmp 
         done 
         cat ${SAMPLE1}_${time_date}/chr*.1M.tmp > ${SAMPLE1}_${time_date}/combine_1M_DP
         #cat ${SAMPLE1}_${time_date}/chr*.1M.DP|cut -f4 > ${SAMPLE1}_${time_date}/combine_1M_DP_un1k
         sh ${path_local}/field_cultivar/script_1M_un1k.sh ${SAMPLE1}_${time_date}
    ) &
wait_all
done 
fi
wait
#-----------------------------------
#-----------------------------------
if true;then
for SAMPLE1 in ${arr[@]};do
      for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
     #1
      ${path_local}/muti_func_snp_compare.py -b 1000000 -i ${SAMPLE1}_${time_date}/${chr}.1M.norm_un1k --mask_cnv on --sample1 ${arr_3[${SAMPLE1}]}  -o ${SAMPLE1}_${time_date}/${chr}.mask_CNV &
     #2
      done
      wait 
      for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
      for item in {deletion,duplication};do
         grep "${item}" ${SAMPLE1}_${time_date}/${chr}.mask_CNV > ${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item} &
      done
      done
wait_all
done
fi
wait
#----------------------------------
#----------------------------------
parallel -j procfile sh bisample_compare_phase1_2.sh ::: $(eval cat sample_path.txt) ::: 2020-04-17162314

#sh bisample_compare_sequence_hugh_sample.sh -t ${time_date} -m ${meta_data} &
         #bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${chr}.2.1M.bed -b /data2/public/ReseqData/workflowFile/${SAMPLE1}/05.mergeAsplit/${chr}.2.*.bam -counts -sorted > ${SAMPLE1}_${time_date}/${chr}.2.1M.DP
