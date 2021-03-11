#!/bin/bash
#set -euxo pipefail
while getopts 'c:v:a:m:t:M:' opt
do
    case "$opt" in
    c) SAMPLE1=$OPTARG ;;
    v) SAMPLE1_name=$OPTARG ;;
    a) SAMPLE1_vcf=$OPTARG ;;
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
if [[ -z "${time_date}" ]]; then
    echo "No time date provided."; exit 1
fi

#SAMPLE_list=${SAMPLE1}
#arr=(`echo ${SAMPLE_list} | tr ',' ' '`)
arr=(`awk '{print $1}' ${SAMPLE1}|tr '\n' ' '|sed s/' '$//g`)
arr_2=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
arr_3=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
path="/data/user/yangzz/mapping/fieldergenomecompare/20200416_201_CNV"
#declare -A arr_2
#while read line;do
#arr_2[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $1}'`
#done < ${meta_data}

#declare -A arr_3
#while read line;do
#arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $4}'`
#done < ${meta_data}

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

for SAMPLE1 in ${arr[@]};do
    #SAMPLE1=`echo $line | awk '{print $1}'`
    for SAMPLE2 in ${arr_2[@]};do
         num=`arritemidx "${arr_2[*]}" ${SAMPLE1}`
         arr_2=(`arrslice "${arr_2[*]}" $num`)
         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
           if [[ ! -d "${SAMPLE1}_${SAMPLE2}" ]]; then
                mkdir ${SAMPLE1}_${SAMPLE2}
           fi 
           (for num in {1..7}; do
             for i in {A,B,D}; do
               for item in {deletion,duplication};do
                 a=`du -h ${path}/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item}`
                 b=`echo $a |awk '{print $1}'`
                 if [ "${b}" == 0 ];then
                   #echo 'file zero' 
                   awk '{print $1"\t"$2"\t"$3"\t"SAMPLE2"'"$item"'""_CNV"}' SAMPLE2="$SAMPLE2" ${path}/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item}  > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item}
                 else
                   #echo 'file not zero' 
                   awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE2"'"$item"'""_CNV"}' SAMPLE2="$SAMPLE2" ${path}/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item} ${path}/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item} > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item}
                 fi
        #
                 a=`du -h ${path}/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item}`
                 b=`echo $a |awk '{print $1}'`
                 #echo ${b}
                 if [ "${b}" == 0 ];then
                    #echo 'file zero'
                    awk '{print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" ${path}/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item}  > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
                 else
                    #echo 'file not zero'
                    awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" ${path}/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item} ${path}/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item} > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
                 fi
                 awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==1) print $1"\t"$2"\t"$3"\t""'"$item"'""_both_CNV"}' ${path}/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item} ${path}/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item} > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}both_CNV_${item}
            #3  
                 cat ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item} ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item} ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}both_CNV_${item} > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV_${item}
               done
               cat ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV_deletion ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV_duplication | sort -u -nk2,2 > ${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV
             done
           done)&
           sleep 0.01s
         fi
#
    done
#wait_all
done
wait
#

#for sample1 in ${arr[@]};do
#for sample2 in ${arr_2[@]};do
#    unset arr_2[${sample1}]
#    if [ "${sample1}" != "${sample2}" ]; then
#    for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do

    #cat tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density_* > tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density
    #./muti_func_snp_compare.py --two_diff_level on -i tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density --sample1 ${sample1} --sample2 ${sample2} --time_data ${time_date} --total ${total} &
#    a1=`du -h tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.${sample1}to${sample2}all_CNV`
#    b1=`echo ${a1} |awk '{print $1}'`
    #cat tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.only_homo_snp_level tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.${sample1}to${sample2}all_CNV > tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_level
    #if [ "${b}" != 0 ];then
    #awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" tmp/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item} tmp/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
    #fi
   #awk '{print $4}' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_density > tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_density
   #sed -i '1i '${sample}_${SAMPLE2}'' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_density
#    done
#    ./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${total}/${sample1}_${sample2}_${time_date}" -s ".homo_snp_level" > plotfile/plotfile_${sample1}_${sample2}
#    Rscript ~/R/graph_diff_sim.R --sample1 ${sample1} --sample1_name ${arr_3[${sample1}]} --sample2 ${sample2} --sample2_name ${arr_3[${sample2}]} -s ${time_date}
#    fi
#paste `awk '{if ($1!="'"${sample}"'")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/""'"${sample}"'""_"$1"_2019-12-15201714/combine_homodiff_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_${sample}_${SAMPLE2}_homodiff_snp_density
#done
#cat `awk '{if ($1!="PH09")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/PH09_"$1"_2019-12-15201714/combine_homo_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_PH09_homo_snp_density
#done




#for sample in ${arr[@]};do
#while read line;do
#    SAMPLE2=`echo $line | awk '{print $1}'`
#    if [ "${sample}" != "${SAMPLE2}" ]; then
#    if [[ -s "tmp/${sample}_${SAMPLE2}_${time_date}/chr1A.homo_snp_level_20" ]]; then
#    cat tmp/${sample}_${SAMPLE2}_${time_date}/chr*.homo_snp_level_* > tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_level
#    awk '{print $5}' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_level > tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_level
#    sed -i '1i '${sample}_${SAMPLE2}'' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_level
#    cp tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_level field_cultivar/tmp_snp_level/${sample}_${SAMPLE2}_homo_snp_level
#    cp tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_level field_cultivar/tmp_snp_level/${sample}_${SAMPLE2}_homodiff_snp_level
#    fi
#    fi
#done < ${meta_data}
#paste `awk '{if ($1!="'"${sample}"'")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/""'"${sample}"'""_"$1"_2019-12-15201714/combine_homodiff_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_${sample}_${SAMPLE2}_homodiff_snp_density
#done
#cat `awk '{if ($1!="PH09")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/PH09_"$1"_2019-12-15201714/combine_homo_snp_density"}' field_cultivar/    metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_PH09_homo_snp_density

