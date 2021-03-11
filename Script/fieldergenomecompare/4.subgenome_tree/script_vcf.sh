#!/bin/bash

# WD=/data2/rawdata2/tree/combineXN_190803/
WD1=/data2/rawdata2/mergeFile/mergeWithXN/A_sub_190807/
WD2=/data2/rawdata2/mergeFile/mergeWithXN/B_sub_190807/
WD3=/data2/rawdata2/mergeFile/mergeWithXN/D_sub_190804

for CHR in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
  bcftools view -M 2 -m 2 -v snps -S /data/user/yangzz/mapping/fieldergenomecompare/1.diff_dev/metadata_cultivar_samplelist.txt --min-ac=1 --threads 4 -a /data/user/yangzz/mapping/08.mergeGVCF/201_final/${CHR}.bcf.gz -Ob -o ${CHR}.trimed.bcf.gz &
  wait_all
done
#
wait
