#!/bin/bash
#set -x
#set -euxo pipefail
while getopts 'v:d:p:m:c:' opt
do
    case "$opt" in
    v) vcf_path=$OPTARG ;;
    d) dev_path=$OPTARG ;;
    p) sample_compare_path=$OPTARG ;;
    m) meta_data=$OPTARG ;;
    c) BED=$OPTARG ;;
    esac
done
if [[ -z "${vcf_path}" ]]; then
    echo "No vcf path provided"; exit 1
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
sample_list=`cut -f2 ${meta_data}|tr '\n' ','|sed  s/','$/''/g`
#dev_path=$2
if false;then
bcftools query -f "[%GT\t]\n"  -H  ${vcf_path}/chr7D.ann.bcf.gz -s  ${sample_list}|head -1 > header.txt
#bcftools query -f "[%GT\t]\n"  -H  /data/user/chenym/100rna-seq/guoyw_leaf_root_h2o2_heterogeneous_system/00.reseq/02.map/total.filter.final.merge.bcf.gz -S  ${sample_list}|head -1 > header.txt
#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do #chr1A
for CHR in `cut -f1 $BED|sort | uniq`; do
n=`cat ${BED} | grep ${CHR} | wc -l`
parallel -j procfile bash ./dist_clus.sh ::: $(eval echo {1..$n}) ::: ${vcf_path}  ::: ${sample_list} ::: ${BED} ::: ${CHR}
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
if false;then
num_use=3
rm -f header_diff.txt
rm -f sample_list.txt
for SAMPLE1 in ${arr[@]};do
    #SAMPLE1=`echo $line | awk '{print $1}'`
    for SAMPLE2 in ${arr_2[@]};do
         num=`arritemidx "${arr_2[*]}" ${SAMPLE1}` 
         arr_2=(`arrslice "${arr_2[*]}" $num`)
         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
         num_use=$(( $num_use + 1 ))
         echo ${SAMPLE1}_${SAMPLE2} >> sample_list.txt
         fi
    done
done
cat sample_list.txt | tr '\n' '\t'| sed s/'\t'$/'\n'/g > header_diff.txt

sed  -i s/'\t'$/'\n'/g chr*.*.diff
mkdir ${dev_path}
mkdir ${dev_path}/diff_dev
mv -f chr*.*.diff ${dev_path}/diff_dev
#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
for CHR in `cut -f1 $BED|sort | uniq`; do
(
cat `sh concat_data.sh ${CHR} ${dev_path}/diff_dev diff $BED` > ${dev_path}/${CHR}.1M.combinediff
cat header_diff.txt ${dev_path}/${CHR}.1M.combinediff |sponge  ${dev_path}/${CHR}.1M.combinediff
) &
wait_all
done
wait
fi
if true;then
#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D  chr3D chr4D chr5D chr6D chr7D;do
for CHR in `cut -f1 $BED|sort | uniq`; do
n=`wc -l < sample_list.txt`
parallel --xapply -j procfile sh cut_file.sh ::: $(eval echo {1..$n}) ::: ${dev_path} ::: ${sample_compare_path} ::: $(eval cat sample_list.txt) ::: ${CHR} ::: ${BED}
done
fi

