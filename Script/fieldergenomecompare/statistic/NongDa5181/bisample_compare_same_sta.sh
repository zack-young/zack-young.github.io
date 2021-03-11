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
if true;then
arr=(`awk '{print $1}' ${SAMPLE1}|tr '\n' ' '|sed s/' '$//g`)
arr_2=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
#arr_3=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
#declare -A arr_2
#while read line;do
#arr_2[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $1}'`#done < ${meta_data}

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
#dev='/data/user/shinyug/HMM_for_yzz_comp/200802/201_sample_compare'
dev='/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare'
for SAMPLE1 in ${arr[@]};do
    #SAMPLE1=`echo $line | awk '{print $1}'`
    for SAMPLE2 in ${arr_2[@]};do
         num=`arritemidx "${arr_2[*]}" ${SAMPLE1}`
         arr_2=(`arrslice "${arr_2[*]}" $num`)
         #echo  "/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare/"${SAMPLE1}_${SAMPLE2} &
         #echo  "/data/user/shinyug/HMM_for_yzz_comp/200802/201_sample_compare/"${SAMPLE1}_${SAMPLE2} &
         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
         if [[ -d "${dev}/${SAMPLE1}_${SAMPLE2}" ]]; then
             sim_num=`cat ${dev}/${SAMPLE1}_${SAMPLE2}/*.homo_hmm_snp_level| grep -E 'low|both'|wc -l`
             awk 'BEGIN{printf "'${SAMPLE1}'""\t""'${SAMPLE2}'""\t""%.2f\n",'${sim_num}'/14075}'>> Similar_region_count_hmm
#           #if [[ ! -d "tmp/${SAMPLE1}_${SAMPLE2}_${time_date}" ]]; then
#           #     mkdir tmp/${SAMPLE1}_${SAMPLE2}_${time_date}
#           #fi 
#         num_tmp=$(( $num_tmp + 1 ))
#         echo $num_tmp
#         for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#         #(./muti_func_snp_compare.py --two_diff_level on -i ${SAMPLE1}_${SAMPLE2}/${chr}.1M.delcnv_density --sample1 ${SAMPLE1} --sample2 ${SAMPLE2} 
          #cat ${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_level ${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV > ${SAMPLE1}_${SAMPLE2}/${chr}.homo_snp_level) &
#         #cat ${SAMPLE1}_${SAMPLE2}/chr*.1M.delcnv_density |grep -v 'CHR' > combine/${SAMPLE1}_${SAMPLE2}_combine.1M.delcnv_density
#         #cat ${SAMPLE1}_${SAMPLE2}/chr*.1M.delcnv_density > combine/header/${SAMPLE1}_${SAMPLE2}_combine_header.1M.delcnv_density
#         #rm ${SAMPLE1}_${SAMPLE2}/combine.1M.delcnv_density
#         #done
#         #( ./ptf_maker.py -p "/data/user/yangzz/mapping/fieldergenomecompare/sample_compare/${SAMPLE1}_${SAMPLE2}" -s ".homo_snp_level" > plotfile/plotfile_${SAMPLE1}_${SAMPLE2}
#         #Rscript ~/R/graph_diff_sim.R --sample1 ${SAMPLE1} --sample1_name ${arr_3[${SAMPLE1}]} --sample2 ${SAMPLE2} --sample2_name ${arr_3[${SAMPLE2}]} -s "2020_2_29" ) &
#         #(cp /data/user/yangzz/mapping/fieldergenomecompare/sample_compare/${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_level > dev/${chr}_${SAMPLE1}_${SAMPLE2}_homo_level) &
#         done
         else
            sim_num=`cat ${dev}/${SAMPLE2}_${SAMPLE1}/*.homo_hmm_snp_level| grep -E 'low|both'|wc -l`
            awk 'BEGIN{printf "'${SAMPLE2}'""\t""'${SAMPLE1}'""\t""%.2f\n",'${sim_num}'/14075}'>> Similar_region_count_hmm
         fi
         fi
    done
#wait_all
done
wait
fi
##
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#cat  dev/${chr}_*_homo_level |grep 'low'> ${chr}_combine_homo_low
#done
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_deletion_split"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_deletion_split
if false; then
thr=3
#pop=`grep -v 'charac' similar_distri/ABD_similar_${thr}_change.csv |awk  '{for(i=2;i<=NF;i++) print $1"_"$i}' |wc -l`
#pop=`wc -l sample_path_hmm.txt`
pop=`wc -l 200_2_nearest_filtered.txt`
for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
        num=${chr:3:2}
        (cat `awk '{print "/data/user/shinyug/HMM_for_yzz_comp/200802/201_sample_compare/"$1"/"chr".only_homo_snp_hmm_level"}' chr="$chr" 200_2_nearest_filtered.txt | tr '\n' ' '`> ${chr}.2_nearest_combine_hmm_level
        #grep 'low' similar_distri/${chr}_similar_up_${thr} |./count_CNV.py --count on --pop_num ${pop} |sort -nk1,1 >  similar_distri/${chr}_similar_up_${thr}_count
        #./count_CNV.py --compensent on --chrom ${num} -i similar_distri/${chr}_similar_up_${thr}_count -o similar_distri/${chr}_similar_up_${thr}_count_compensent
        echo $num
        grep 'low' ${chr}.2_nearest_combine_hmm_level |./count_CNV.py --count on --pop_num ${pop} |sort -nk1,1 >  similar_distri/${chr}.2_nearest_combine_hmm_level_count
        ./count_CNV.py --compensent on --chrom ${num} -i similar_distri/${chr}.2_nearest_combine_hmm_level_count -o similar_distri/${chr}.2_nearest_combine_hmm_level_count_compensent) &
done
fi
