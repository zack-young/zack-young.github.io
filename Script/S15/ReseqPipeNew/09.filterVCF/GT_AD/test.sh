#!/bin/bash
while getopts 'c:m:' opt
do
    case "$opt" in
    c) SAMPLE1=$OPTARG ;;
    m) meta_data=$OPTARG ;;
    esac
done
if [[ -z "${SAMPLE1}" ]]; then
    echo "No first sample provided"; exit 1
fi
if [[ -z "${meta_data}" ]]; then
    echo "No sample data provided."; exit 1
fi

a=${SAMPLE1}
arr=(`echo ${a} | tr ',' ' '`)
echo ${arr[@]}
while read line;do
SAMPLE1=`echo $line | awk '{print $1}'`
num=0
echo ${SAMPLE1}
for SAMPLE2 in ${arr[@]};do
    if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
    b=1
    #echo ${SAMPLE1} ${SAMPLE2}
    else
    echo ${SAMPLE2}
    echo ${arr[${num}]}
    echo ${num}
    unset arr[${num}]
    arr=(`echo ${arr[@]}`)
    echo ${arr[@]}
    fi
    num=$(($num + 1))
done
done < ${meta_data}

    
