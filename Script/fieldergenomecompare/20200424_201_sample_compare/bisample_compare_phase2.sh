#!/bin/bash
while getopts 'c:m:t:p:' opt
do
    case "$opt" in
    c) SAMPLE1=$OPTARG ;;
    m) meta_data=$OPTARG ;;
    t) time_date=$OPTARG ;;
    p) DEV_PATH=$OPTARG ;;
    esac
done
if [[ -z "${SAMPLE1}" ]]; then
    echo "No first sample provided"; exit 1
fi
if [[ -z "${meta_data}" ]]; then
    echo "No sample data provided."; exit 1
fi
if [[ -z "${DEV_PATH}" ]]; then
    echo "No CNV file path provided."; exit 1
fi

#SAMPLE_list=${SAMPLE1}
#arr=(`echo ${SAMPLE_list} | tr ',' ' '`)
arr=(`awk '{print $1}' ${SAMPLE1}|tr '\n' ' '|sed s/' '$//g`)
arr_2=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
arr_3=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
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

for SAMPLE1 in ${arr[@]};do
    #SAMPLE1=`echo $line | awk '{print $1}'`
    for SAMPLE2 in ${arr_2[@]};do
         num=`arritemidx "${arr_2[*]}" ${SAMPLE1}`
         arr_2=(`arrslice "${arr_2[*]}" $num`)
         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
           #if [[ ! -d "tmp/${SAMPLE1}_${SAMPLE2}_${time_date}" ]]; then
           #     mkdir tmp/${SAMPLE1}_${SAMPLE2}_${time_date}
           #fi 
           (for num in {1..7}; do
             for i in {A,B,D}; do
                 a=`du -h ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV`
                 b=`echo $a |awk '{print $1}'`
                 if [ "${b}" == 0 ];then
                   awk '{print $0}'   ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.density  > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.delcnv_density
                   #awk '{print $0}'   ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.hete_miss_density  > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.hete_miss_delcnv_density
                 else
                   awk 'NR==FNR{a[$2]=$2;next}NR>FNR{if($2 in a ==0) print $0}'  ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV  ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.density  > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.delcnv_density
                   #awk 'NR==FNR{a[$2]=$2;next}NR>FNR{if($2 in a ==0) print $0}'  ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV  ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.hete_miss_density  > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.1M.hete_miss_delcnv_density
                 fi
        #
             done
           done) &
           #wait
           #cat ${SAMPLE1}_${SAMPLE2}/chr*.1M.delcnv_density > statistic/combine.${SAMPLE1}_${SAMPLE2}_1M.delcnv_density
         fi
         sleep 0.01s
#
    done
wait_all
done
wait
#
