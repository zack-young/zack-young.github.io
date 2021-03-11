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
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $2}'`
done < ${meta_data}
#echo ${arr_3[@]}
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
    num=0
    #echo "/data/user/yangzz/mapping/fieldergenomecompare/statistic/1B1R_lovrin/chr1B.${SAMPLE1}_undefined_level	chr1B	${arr_3[${SAMPLE1}]}"
    #SAMPLE1=`echo $line | awk '{print $1}'`
    for SAMPLE2 in SRR10766512 S212;do  #${arr_2[@]};do
         #num=`arritemidx "${arr_2[*]}" ${SAMPLE1}`
         #arr_2=(`arrslice "${arr_2[*]}" $num`)
         #SAMPLE2='S212'
         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
         #echo  "/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare/"${SAMPLE1}_${SAMPLE2} &
         #echo  "/data/user/shinyug/HMM_for_yzz_comp/200802/201_sample_compare/"${SAMPLE1}_${SAMPLE2} &
         if [[ -d "${dev}/${SAMPLE1}_${SAMPLE2}" ]]; then
         echo "${dev}/${SAMPLE1}_${SAMPLE2}/chr1B.snp_level_sp	NA	${num}"         
         #sort -k 2,2 -n ${dev}/${SAMPLE1}_${SAMPLE2}/chr1B.homo_undefined_snp_level | sponge ${dev}/${SAMPLE1}_${SAMPLE2}/chr1B.homo_undefined_snp_level
         #sed s/'low'/"low${num}"/g ${dev}/${SAMPLE1}_${SAMPLE2}/chr1B.homo_undefined_snp_level > ${dev}/${SAMPLE1}_${SAMPLE2}/chr1B.snp_level_sp
         #array_name[$num]=`echo  "/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare/"${SAMPLE1}_${SAMPLE2}`
         else
         echo "${dev}/${SAMPLE2}_${SAMPLE1}/chr1B.snp_level_sp	NA	${num}"
         #sort -k 2,2 -n ${dev}/${SAMPLE2}_${SAMPLE1}/chr1B.homo_undefined_snp_level | sponge ${dev}/${SAMPLE2}_${SAMPLE1}/chr1B.homo_undefined_snp_level
         #sed s/'low'/"low${num}"/g ${dev}/${SAMPLE2}_${SAMPLE1}/chr1B.homo_undefined_snp_level > ${dev}/${SAMPLE2}_${SAMPLE1}/chr1B.snp_level_sp
         #array_name[$num]=`echo  "/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare/"${SAMPLE2}_${SAMPLE1}`
         fi
         num=$(( $num + 1))
#         cat ${dev}/${SAMPLE1}_S212/chr1B.homo_hmm_snp_level |sort -nk 2,2 > chr1B.${SAMPLE1}_homo_hmm_snp_level
#         else
#         cat ${dev}/S212_${SAMPLE1}/chr1B.homo_hmm_snp_level| sort -nk 2,2 >  chr1B.${SAMPLE1}_homo_hmm_snp_level
         #cut -f4 chr1B.${SAMPLE1}_homo_hmm_snp_level | paste chr1B.lovrin_simmilar.dist - | sponge chr1B.lovrin_simmilar.dist
         #same_num=`grep -E 'low|both' ${dev}/${SAMPLE1}_${SAMPLE2}/chr1B.homo_hmm_snp_level|wc -l`
         #all_same_num=`cat ${dev}/${SAMPLE1}_${SAMPLE2}/chr*.homo_hmm_snp_level | grep -E 'low|both'|wc -l`
         #echo "${same_num} ${all_same_num} ${arr_3[${SAMPLE1}]} ${SAMPLE1}" >> sim_count
         #elif [[ -d "${dev}/${SAMPLE2}_${SAMPLE1}" ]]; then
         #same_num=`grep -E 'low|both' ${dev}/${SAMPLE2}_${SAMPLE1}/chr1B.homo_hmm_snp_level|wc -l`
         #all_same_num=`cat ${dev}/${SAMPLE2}_${SAMPLE1}/chr*.homo_hmm_snp_level | grep -E 'low|both'|wc -l`
         #echo "${same_num} ${all_same_num} ${arr_3[${SAMPLE1}]} ${SAMPLE1}" >> sim_count
         #fi
#         num_tmp=$(( $num_tmp + 1 ))
#         echo $num_tmp
         fi
    done
echo "NA	${arr_3[${SAMPLE1}]}	NA"
#echo ${array_name[*]}
#paste ${array_name[0]}/chr1B.homo_undefined_snp_level  ${array_name[1]}/chr1B.homo_undefined_snp_level| cut -f1,2,3,4,8|awk '{if ($0!~"low"&& $0!~"both" ){print $1"\t"$2"\t"$3"\t""undefined"}}' > chr1B.${SAMPLE1}_undefined_level
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
