#!/bin/bash
#set -euxo pipefail
while getopts 'c:v:a:d:m:p:' opt
do
    case "$opt" in
    c) SAMPLE1=$OPTARG ;;
    v) SAMPLE1_name=$OPTARG ;;
    a) SAMPLE1_vcf=$OPTARG ;;
    d) total=$OPTARG ;;
    m) meta_data=$OPTARG ;;
    p) DEV_PATH=$OPTARG ;;
    esac
done
if [[ -z "${SAMPLE1}" ]]; then
    echo "No first sample provided"; exit 1
fi
if [[ -z "${meta_data}" ]]; then
    echo "No sample data provided."; exit 1
fi
if false;then
#SAMPLE_list=${SAMPLE1}
#arr=(`echo ${SAMPLE_list} | tr ',' ' '`)
arr=(`awk '{print $1}' ${SAMPLE1}|tr '\n' ' '|sed s/' '$//g`)
arr_2=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
#arr_3=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
#declare -A arr_2
#while read line;do
#arr_2[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $1}'`
#done < ${meta_data}

declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $4}'`
done < ${meta_data}

#for sample1 in ${arr[@]};do
#for sample2 in ${arr_2[@]};do
#    unset arr_2[${sample1}]
#    if [ "${sample1}" != "${sample2}" ]; then
arritemidx(){
  local tmp
  local count=0
  local array=`echo $1`
  for tmp in ${array[@]};do
    if test $2 = $tmp;then
      echo $count
      return
    fi
    count=$(( $count + 1 ))
  done
  echo -1
}

arrslice(){
  array=($1)
  if [ $2 == -1 ];then
    echo ${array[@]}
  elif [ $2 == 0 ];then
    echo ${array[@]:1}
  else
    #echo "${array[@]:0:$2} ${array[@]:$(( $2 + 1 ))}"
    echo "${array[*]:0:$2} ${array[*]:$(( $2 + 1 ))}"
  fi
}
num_tmp=0
line_num=1
#for SAMPLE1 in ${arr[@]};do
#    for SAMPLE2 in ${arr_2[@]};do
#         num=`arritemidx "${arr_2[*]}" ${SAMPLE1}`
#         arr_2=(`arrslice "${arr_2[*]}" $num`)
fi
while read line;do
    SAMPLE1=`echo ${line}|awk '{print $1}'`
    SAMPLE2=`echo ${line}|awk '{print $2}'`
         (for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
             if [[ -d "${DEV_PATH}/${SAMPLE1}_${SAMPLE2}" ]]; then
#             if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
             num_tmp=$(( $num_tmp + 1 ))
             #echo $num_tmp
             #echo ${SAMPLE1}_${SAMPLE2}
             #cat ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/${chr}.homo_hmm_snp_level >> ${chr}.4_sim_combine_hmm_level
             ./count_CNV.py --crossover on -i ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/${chr}.homo_hmm_snp_level >> ${chr}.homo_snp_hmm_crossover
             else
             #cat ${DEV_PATH}/${SAMPLE2}_${SAMPLE1}/${chr}.homo_hmm_snp_level >> ${chr}.4_sim_combine_hmm_level
             #echo ${SAMPLE2}_${SAMPLE1}
             ./count_CNV.py --crossover on -i ${DEV_PATH}/${SAMPLE2}_${SAMPLE1}/${chr}.homo_hmm_snp_level >> ${chr}.homo_snp_hmm_crossover
             fi
         done) &
         sleep 0.01s
#    done
#echo $line_num
#line_num=$(( $line_num + 1 ))
#done
done < ${meta_data} 
#fi

##
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#cat  dev/${chr}_*_homo_level |grep 'low'> ${chr}_combine_homo_low
#done
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_deletion_split"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_deletion_split
if false;then
pop=`wc -l 200_2_nearest_filtered.txt`
for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
    num=${chr:3:2}
    (#./count_CNV.py --crossover on -i ${chr}.2_nearest_combine_hmm_level -o ${chr}.2_nearest_homo_snp_hmm_crossover
    ./count_CNV.py --count on --pop_num ${pop} -i ${chr}.2_nearest_homo_snp_hmm_crossover|sort -nk1,1 >  crossover_distri/${chr}.2_nearest_combine_hmm_crossover_count
    ./count_CNV.py --compensent on --chrom ${num} -i crossover_distri/${chr}.2_nearest_combine_hmm_crossover_count -o crossover_distri/${chr}.2_nearest_combine_hmm_crossover_count_compensent
    ) &
        #./count_CNV.py --count on --pop_num 7021 -i ${chr}_combine_homo_low |sort -nk1,1 >  ${chr}_combine_homo_low_count
        #./count_CNV.py --compensent on --chrom ${num} -i ${chr}_combine_homo_low_count -o ${chr}.combine_homo_low_compensent
done
fi
