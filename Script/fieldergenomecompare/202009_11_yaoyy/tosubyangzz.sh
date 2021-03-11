#!/bin/bash


#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#  n=$(wc -l < ${CHR}.1M.bed)
#  parallel -j procfile bash ./dist_clus.sh ${CHR} {} ::: $(eval echo {1..$n})
#done

#sed  -i s/'\t'$/'\n'/g chr*.*.diff
#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#cat `sh concat_data.sh ${CHR}` > ${CHR}.1M.combinediff
#done


meta_data=$1
SAMPLE_list=`awk '{print $1}' ${meta_data}|tr '\n' ','|sed s/,$//g`
arr=(`echo ${SAMPLE_list} | tr ',' ' '`)
arr2=(`echo ${SAMPLE_list} | tr ',' ' '`)

#for i in "${!arr[@]}";do
#for a in "${!arr2[@]}";do
#unset arr2[${i}]
#if [ "${i}" != "${a}" ]; then
#echo "${arr[$i]}_${arr2[$a]}"
#fi
#done
#done
#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#paste -d '\t' ${CHR}.1M.bed.tmp ${CHR}.1M.combinediff.tmp > ${CHR}.1M.combinediff
#done

CHR=$2

num=3
for i in "${!arr[@]}";do
for a in "${!arr2[@]}";do
unset arr2[${i}]
if [ "${i}" != "${a}" ]; then
if [[ ! -d "${arr[$i]}_${arr2[$a]}" ]]; then
     mkdir ${arr[$i]}_${arr2[$a]}
fi
num=$[$num+1]
cut -f${num} ../dev/${CHR}.1M.combinediff > ${arr[$i]}_${arr2[$a]}/${CHR}.1M.density 
paste -d '\t' ../dev/${CHR}.1M.bed.tmp ${arr[$i]}_${arr2[$a]}/${CHR}.1M.density | sponge ${arr[$i]}_${arr2[$a]}/${CHR}.1M.density
fi
done
done
