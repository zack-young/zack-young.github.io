#!/usr/bin/env bash
#set -x

#while getopts 'g:a:f:n:' opt
#do
#    case "$opt" in
#    g) GENE=$OPTARG ;;
#    a) LABEL=$OPTARG ;;
#    f) FLANK=$OPTARG ;;
#    n) NONE=$OPTARG ;;
#    esac
#done
#
#if [[ -z "${GENE}" ]]; then
#    echo "No gene provided."; exit 1
#fi
#if [[ ! -e "${LABEL}" ]]; then
#    echo "No annotation file provided"; exit 2
#fi
#if [[ -z "${FLANK}" ]]; then
#    FLANK=0
#fi
#if [[ ! -e "${NONE}" ]]; then
#    NONE="none-var-gene.txt"
#fi
#
#shift $[($OPTIND-1)]
#echo "there are $# group files"
#
#for ((i=1; i<=$#; i++)); do
#    #printf "%d\t%s\n" $i "${!i}"
#    Group=$(basename "${!i}"|cut -d. -f1)
#    gawk -vG=${Group} '{print $1"\t"G}' "${!i}" 
#    gawk -vG=${Group} '{print $1}' "${!i}" 
#done
#arr=(23 56 99 22 11)
#unset arr[0]
#arr1=(`echo ${arr[@]}`)
#echo ${arr[@]}
#echo ${arr1[1]}
declare -A arr2
while read line;do
arr2[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $4}'`
done < field_cultivar/metadata_closersample.txt
sample='S142'
echo ${arr2[${sample}]}

#for sample in ${arr2[@]};do
#for sample2 in ${arr2[@]};do
#unset arr2[${sample}]
#if [ "${sample}" != "${sample2}" ]; then
#echo "${sample},${sample2}"
#fi
#done
#done


