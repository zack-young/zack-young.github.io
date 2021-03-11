#!/bin/bash

# WD=/data2/rawdata2/tree/combineXN_190803/
WD1=/data2/rawdata2/mergeFile/mergeWithXN/A_sub_190807/
WD2=/data2/rawdata2/mergeFile/mergeWithXN/B_sub_190807/
WD3=/data2/rawdata2/mergeFile/mergeWithXN/D_sub_190804

for CHR in chr1A chr2A chr3A chr4A chr5A chr6A chr7A;do
  bcftools view -M 2 -m 2 -v snps -S AABBDDsample_list.txt --min-ac=1 --threads 4 -a ${WD1}/${CHR}.bcf.gz -Ob -o ${CHR}.trimed.bcf.gz &
  wait_all
done
#
for CHR in chr1B chr2B chr3B chr4B chr5B chr6B chr7B;do
  bcftools view -M 2 -m 2 -v snps -S AABBDDsample_list.txt --min-ac=1 --threads 4 -a ${WD2}/${CHR}.bcf.gz -Ob -o ${CHR}.trimed.bcf.gz &
  wait_all
done
#
for CHR in chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
  bcftools view -M 2 -m 2 -v snps -S AABBDDsample_list.txt --min-ac=1 --threads 4 -a ${WD3}/${CHR}.bcf.gz -Ob -o ${CHR}.trimed.bcf.gz &
  wait_all
done

wait
