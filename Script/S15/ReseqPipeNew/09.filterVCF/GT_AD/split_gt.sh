#!/bin/bash
set -x

for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
#for chr in chr1A chr1B chr1D;do
pos=`tail -1 ${chr}.snp_matrix | awk '{print $2}' | awk '{print int($1/1000000)+1}' `
line=`tail -1 ${chr}.snp_matrix | awk '{print $2}' | awk '{print int($1/20000000)+1}' `

#for i in `seq 9 10 | tr '\n' ' '|sed s/,$//g`;do
for i in `seq 1 20 | tr '\n' ' '|sed s/,$//g`;do
left=$[0+${line}*(${i}-1)]
right=$[${left}+$line]
echo $left $right
cat vcf_gt_header > ${chr}.snp_matrix_${i}
awk 'left*1000000 <= $2 && $2 <= right*1000000 {print $0}' left="$left" right="$right"  ${chr}.snp_matrix >>  ${chr}.snp_matrix_${i} &
done
wait_all
done
