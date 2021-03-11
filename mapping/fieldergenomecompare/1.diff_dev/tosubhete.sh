#!/bin/bash
#set -x
#set -euxo pipefail
while getopts 'v:s:d:p:m:' opt
do
    case "$opt" in
    v) vcf_path=$OPTARG ;;
    s) sample_list=$OPTARG ;;
    d) dev_path=$OPTARG ;;
    p) sample_compare_path=$OPTARG ;;
    m) meta_data=$OPTARG ;;
    esac
done
if [[ -z "${vcf_path}" ]]; then
    echo "No vcf path provided"; exit 1
fi
if [[ -z "${sample_list}" ]]; then
    echo "No sample list provided."; exit 1
fi
if [[ -z "${dev_path}" ]]; then  #store diff file
    echo "No dev path provided."; exit 1
fi
if [[ -z "${sample_compare_path}" ]]; then #store sample cross compare file
    echo "No sample compare path provided."; exit 1
fi
if [[ -z "${meta_data}" ]]; then
    echo "No meta data provided."; exit 1
fi

#vcf_path="/data/user/yangzz/mapping/08.mergeGVCF/field_cultivar"
#sample_list="metadata_cultivar_samplelist.txt"
#dev_path=$2
if true;then
bcftools query -f "[%GT\t]\n"  -H  ${vcf_path}/chr1A.bcf.gz -S  ${sample_list}|head -1 > header_hete.txt
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#for CHR in chr3A; do
  n=$(wc -l < ${CHR}.1M.bed)
  parallel -j procfile bash ./dist_clus_hete.sh ${CHR} {} ::: $(eval echo {1..$n}) ::: ${vcf_path}  ::: ${sample_list}
done
fi

#meta_data=$1
SAMPLE_list=`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`
arr=(`echo ${SAMPLE_list}`)
arr_2=(`echo ${SAMPLE_list}`)

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
if true;then
num_use=3
rm header_diff_hete.txt
for SAMPLE1 in ${arr[@]};do
    #SAMPLE1=`echo $line | awk '{print $1}'`
    for SAMPLE2 in ${arr_2[@]};do
         num=`arritemidx "${arr_2[*]}" ${SAMPLE1}` 
         arr_2=(`arrslice "${arr_2[*]}" $num`)
         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
         num_use=$(( $num_use + 1 ))
         echo ${SAMPLE1}_${SAMPLE2} >> header_diff_hete.txt
         fi
    done
done
cat header_diff_hete.txt | tr '\n' '\t'| sed s/'\t'$/'\n'/g | sponge header_diff_hete.txt

#sed  -i s/'\t'$/'\n'/g chr*.*.only_miss_diff
sed  -i s/'\t'$/'\n'/g chr*.*.hete_miss_diff
mkdir ${dev_path}/hete_miss_dev
#mkdir ${dev_path}/only_miss_dev
mv chr*.*.hete_miss_diff ${dev_path}/hete_miss_dev
#mv chr*.*.only_miss_diff ${dev_path}/only_miss_dev
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#for CHR in chr3A; do
(
cat `sh concat_data.sh ${CHR} ${dev_path}/hete_miss_dev hete_miss_diff` > ${dev_path}/${CHR}.1M.combine_hete_miss_diff
cat header_diff_hete.txt ${dev_path}/${CHR}.1M.combine_hete_miss_diff |sponge  ${dev_path}/${CHR}.1M.combine_hete_miss_diff
#cat `sh concat_data.sh ${CHR} ${dev_path}/only_miss_dev only_miss_diff` > ${dev_path}/${CHR}.1M.combine_only_miss_diff
#cat header_diff_hete.txt ${dev_path}/${CHR}.1M.combine_only_miss_diff |sponge  ${dev_path}/${CHR}.1M.combine_only_miss_diff
) &
wait_all
done
fi
if false;then
arr=(`echo ${SAMPLE_list}`)
arr_2=(`echo ${SAMPLE_list}`)
num_use=0
for SAMPLE1 in ${arr[@]};do
    #SAMPLE1=`echo $line | awk '{print $1}'`
    for SAMPLE2 in ${arr_2[@]};do
         num=`arritemidx "${arr_2[*]}" ${SAMPLE1}`
         arr_2=(`arrslice "${arr_2[*]}" $num`)
         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then

         if [[ ! -d "${sample_compare_path}/${SAMPLE1}_${SAMPLE2}" ]]; then
           mkdir ${sample_compare_path}/${SAMPLE1}_${SAMPLE2}
         fi
         num_use=$(( $num_use + 1 ))
         echo ${num_use}
         for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do 
         #for CHR in chr3A; do
         (
         cut -f${num_use} ${dev_path}/${CHR}.1M.combine_hete_miss_diff > ${sample_compare_path}/${SAMPLE1}_${SAMPLE2}/${CHR}.1M.hete_miss_density
          
         paste -d '\t' ${CHR}.1M.bed.tmp ${sample_compare_path}/${SAMPLE1}_${SAMPLE2}/${CHR}.1M.hete_miss_density | sponge ${sample_compare_path}/${SAMPLE1}_${SAMPLE2}/${CHR}.1M.hete_miss_density
         #cut -f${num_use} ${dev_path}/${CHR}.1M.combine_diff > ${sample_compare_path}/${SAMPLE1}_${SAMPLE2}/${CHR}.1M.density
         #paste -d '\t' ${CHR}.1M.bed.tmp ${sample_compare_path}/${SAMPLE1}_${SAMPLE2}/${CHR}.1M.density | sponge ${sample_compare_path}/${SAMPLE1}_${SAMPLE2}/${CHR}.1M.density
         ) &
         done
         fi
    done
    wait_all
done
fi
