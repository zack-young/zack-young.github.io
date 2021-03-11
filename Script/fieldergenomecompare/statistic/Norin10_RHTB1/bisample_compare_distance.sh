#!/bin/bash
#set -euxo pipefail
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
if [[ -z "${meta_data}" ]]; then
    echo "No sample data provided."; exit 1
fi

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

for SAMPLE1 in ${arr[@]};do
    #SAMPLE1=`echo $line | awk '{print $1}'`
    #for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
    chr=chr4B
    num_low=0
    num_diff=0
    for num in `seq 30000001 1000000 31000001`
      if [[ ! -s "../homo_level_dev/${chr}.S238_${SAMPLE1}_homo_snp_level" ]]; then
      item=`awk '{if($2==num) print $4}' num="$num" ${chr}.S238_${SAMPLE1}_homo_snp_level`
      else
      item=`awk '{if($2==num) print $4}' num="$num" ${chr}.${SAMPLE1}_S238_homo_snp_level`
      fi
      if [ "$item" = "low" ];then
        #num_tmp=$(($num_tmp + 1))
        echo -e "${SAMPLE1}\t1" >> ${chr}_sim.txt
      elif [ "$item" = "mid" ]||[ "$item" = "high" ] ;then
        #num_diff=$(($num_diff + 1))
        echo -e "${SAMPLE1}\t0" >> ${chr}_sim.txt
      fi
      #final_num=`expr $num_tmp/($num_tmp+$num_diff)`
    chr=chr4D
    num=18000001
    if [[ ! -s "../homo_level_dev/${chr}.S238_${SAMPLE1}_homo_snp_level" ]]; then
    item=`awk '{if($2==num) print $4}' num="$num" ${chr}.S238_${SAMPLE1}_homo_snp_level`
    else
    item=`awk '{if($2==num) print $4}' num="$num" ${chr}.${SAMPLE1}_S238_homo_snp_level`
    fi
    if [ "$item" = "low" ];then
       echo -e "${SAMPLE1}\t2" >> ${chr}_sim.txt
    else
       echo -e "${SAMPLE1}\t0" >> ${chr}_sim.txt
    fi
#wait_all
done
wait
cut -f2 chr4D_sim.txt|paste chr4B_sim.txt -|awk '{print $1"\t"$2+$3 }'  >chr4B_4D_sim.txt
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

