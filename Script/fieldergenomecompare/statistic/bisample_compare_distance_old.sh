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
for SAMPLE1 in ${arr[@]};do
    #SAMPLE1=`echo $line | awk '{print $1}'`
    #
    #for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
    #    echo ${SAMPLE1} > distance_matrix/${chr}_${SAMPLE1}_torestdis
    #    for((i=1;i<=${line_num};i++));do   
    #        echo -en '\n' >> distance_matrix/${chr}_${SAMPLE1}_torestdis
    #    done
    #done
    for SAMPLE2 in ${arr_2[@]};do
    #     num=`arritemidx "${arr_2[*]}" ${SAMPLE1}`
    #     arr_2=(`arrslice "${arr_2[*]}" $num`)
         (for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
             #if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
             if [[ -d "${DEV_PATH}/${SAMPLE1}_${SAMPLE2}" ]]; then
             num_tmp=$(( $num_tmp + 1 ))
             #echo $num_tmp
             numerator1=`grep 'low' ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/${chr}.homo_hmm_snp_level|wc -l`
             numerator2=`grep 'both_CNV' ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/${chr}.homo_hmm_snp_level|wc -l`
             denominator=`cat ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/${chr}.homo_hmm_snp_level | wc -l `
             decimal=`awk 'BEGIN{printf "%.4f\n",('$numerator1'+'$numerator2')/'$denominator'}'`
             echo $decimal > distance_matrix/${chr}_${SAMPLE1}_${SAMPLE2}_torestdis
             #(cp /data/user/yangzz/mapping/fieldergenomecompare/sample_compare/${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_level > dev/${chr}_${SAMPLE1}_${SAMPLE2}_homo_level) &
             #done
             else
             echo 0 > distance_matrix/${chr}_${SAMPLE1}_${SAMPLE2}_torestdis
             fi
         done) &
         sleep 0.01s
    done
echo $line_num
line_num=$(( $line_num + 1 ))
done
wait

for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#for chr in chr6A chr6B;do
while read line;do
sample=`echo $line | awk '{print $1}'`
cat `awk '{print "distance_matrix/"chr"_"sample"_"$1"_torestdis"}' chr="$chr" sample="$sample" ${meta_data}|tr '\n' ' '|sed s/,$//g` > distance_matrix/chr_${sample}_torestdis
done < ${meta_data}
paste `awk '{print "distance_matrix/"chr"_"$1"_torestdis"}' chr="$chr" ${meta_data} |tr '\n' ' '|sed s/,$//g` > ${chr}_200_distance
awk '{print $1}'  ${meta_data}|tr '\n' '\t'|sed s/"\t"$/"\n"/g > header.txt
cat header.txt ${chr}_200_distance | sponge ${chr}_200_distance
#cat `grep -v 'charac'  ABD_similar_2_change.csv | awk  '{for(i=2;i<=NF;i++) print "dev/"chr"_"$1"_"$i"_homo_level"}' chr="$chr" |tr '\n' ' '|sed s/,$//g` > ${chr}_similar_up_2 &
#cat `grep -v 'charac'  ABD_similar_3_change.csv | awk  '{for(i=2;i<=NF;i++) print "dev/"chr"_"$1"_"$i"_homo_level"}' chr="$chr" |tr '\n' ' '|sed s/,$//g` > ${chr}_similar_up_3 &
#cp `grep -v 'charac' awk 'BEGIN{printf "%.2f\n",('$numerator'+'$numerator')/'$denominator'}'  similar_distri/ABD_similar_1_change.csv | awk  '{for(i=2;i<=NF;i++) print "dev/"chr"_"$1"_"$i"_homo_level"}' chr="$chr" `  similar_distri/length/

done

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

