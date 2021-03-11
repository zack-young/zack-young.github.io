#!/bin/bash
set -euxo pipefail

while getopts 'l:b:' opt
do
    case "$opt" in
    l) meta_data=$OPTARG ;;
    b) BED=$OPTARG ;;
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


if false;then
for SAMPLE1 in ${arr[@]};do
    (if [[ ! -d "${SAMPLE1}_${time_date}" ]]; then
         mkdir ${SAMPLE1}_${time_date}
     fi
         for CHR in `cut  -f1 ${BED}|sort | uniq`; do
         grep ${CHR} ${BED} > ${CHR}_${SAMPLE1}.bed
         awk '{if ($2>400000000) {print CHR".2""\t"$2-400000000"\t"$3-400000000} else {print CHR".1""\t"$2"\t"$3}}' CHR="$CHR" ${CHR}_${SAMPLE1}.bed > ${CHR}_${SAMPLE1}.bed2
         bedtools coverage -a ${CHR}_${SAMPLE1}.bed2 -b /data/user/yangzz/mapping/Reseq_data/${SAMPLE1}/05.mergeAsplit/${CHR}*.bam -counts -sorted > ${SAMPLE1}_${time_date}/${CHR}.1M.DP
         gawk '{print $4*1000/($3-$2)}' ${SAMPLE1}_${time_date}/${CHR}.1M.DP > ${SAMPLE1}_${time_date}/${CHR}.1M.tmp
         sh script_1M_un1k.sh dev/${SAMPLE1}_${time_date}  ${SAMPLE1}_${time_date} ${CHR}
         cut -f1,2,3 ${CHR}_${SAMPLE1}.bed | paste - ${SAMPLE1}_${time_date}/${CHR}.1M.norm_un1k | sponge ${SAMPLE1}_${time_date}/${CHR}.1M.norm_un1k
         rm -f ${CHR}_${SAMPLE1}.bed
         rm -f ${CHR}_${SAMPLE1}.bed2
         done
    ) &
wait_all
done 
fi
wait
#-----------------------------------
#-----------------------------------
if false;then
for SAMPLE1 in ${arr[@]};do
      for chr in `cut -f1 $BED|sort | uniq`; do
     #1
      ./muti_func_snp_compare.py -b 1000000 -i ${SAMPLE1}_${time_date}/${chr}.1M.norm_un1k --mask_cnv on --sample1 ${arr_3[${SAMPLE1}]}  -o ${SAMPLE1}_${time_date}/${chr}.mask_CNV &
     #2
      done
      wait 
      for chr in `cut -f1 $BED|sort | uniq`; do
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
parallel -j procfile sh bisample_compare_phase1_2.sh ::: $(eval cat sample_path.txt) ::: 2020-04-17162314 ::: ${BED}

#sh bisample_compare_sequence_hugh_sample.sh -t ${time_date} -m ${meta_data} &
         #bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${chr}.2.1M.bed -b /data2/public/ReseqData/workflowFile/${SAMPLE1}/05.mergeAsplit/${chr}.2.*.bam -counts -sorted > ${SAMPLE1}_${time_date}/${chr}.2.1M.DP
