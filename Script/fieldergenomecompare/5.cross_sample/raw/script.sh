#!/bin/bash

CHR=$1
N=$2
WD="/data/user/yangzz/mapping/08.mergeGVCF/201_final"
PA="/data/user/yangzz/mapping/fieldergenomecompare/1.diff_dev/metadata_cultivar_samplelist.txt"
/data/user/yangzz/worktools/Plink/plink --bcf <(bcftools view -v snps --min-ac=1 -M2 -m2 -R <(tail -n +${N} ${CHR}.1M.bed|head -n 1) ${WD}/${CHR}.ann.bcf.gz -S ${PA} -Ob) \
      --allow-extra-chr \
      --distance square flat-missing \
      --out ${CHR}.${N}.rawdist
