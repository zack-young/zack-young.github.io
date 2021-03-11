#!/bin/bash
#set -euxo pipefail
meta_data=$1
SAMPLE1=$2
dev_path=$3

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
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
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
#dir="/data/user/yangzz/mapping/fieldergenomecompare/statistic/Rht8"
if true;then
for SAMPLE1 in ${arr[@]};do
    for SAMPLE2 in ${arr_2[@]};do
         num=`arritemidx "${arr_2[*]}" ${SAMPLE1}`
         arr_2=(`arrslice "${arr_2[*]}" $num`)
         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
           num_tmp=$(( $num_tmp + 1 ))
           num_low=0
           num_diff=0
           num_cnv=10
           #echo $num_tmp
           #if [[ -d "${dev_path}/${SAMPLE1}_${SAMPLE2}" ]]; then
                                                             
           echo ${SAMPLE1},${SAMPLE2}
           #else
           #  echo ${SAMPLE2},${SAMPLE1}
           #fi
         fi
    done
#wait_all
line_num=$(( $line_num + 1 ))
done
fi
#paste `awk '{print "'$dir'""/""'$chr'""_"$1"_""'$suffix2'""_dist"}'  ../metadata_cultivar_final.txt |tr '\n' ' '|sed s/,$//g` > ${dir}/${chr}_combine_${suffix2}_dist
##
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#cat  dev/${chr}_*_homo_level |grep 'low'> ${chr}_combine_homo_low
#done
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_deletion_split"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_deletion_split
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
        #num=${chr:3:2}
        #./count_CNV.py --count on --pop_num 7021 -i ${chr}_combine_homo_low |sort -nk1,1 >  ${chr}_combine_homo_low_count
        #./count_CNV.py --compensent on --chrom ${num} -i ${chr}_combine_homo_low_count -o ${chr}.combine_homo_low_compensent
#done

