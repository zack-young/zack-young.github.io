#!/usr/bin/env bash
set -euxo pipefail

PLOTSCRIPT=/data2/rawdata2/tree/script/NJtree_iTOL.R
COLORFILE=/data2/rawdata2/tree/combineXN_190803/colorF_5cols.txt
SAMPLEORDER=/data2/rawdata2/tree/combineXN_190803/D_sub_190809/03.plot/sample_order_hex.txt

for maf in 01 02 03 04 05 1;do
for miss in 05 1 2 3 4 5;do
(for CHR in AABB DD AABBDD;do
  (#plink --bfile /data2/rawdata2/tree/combineXN_190803/rawdata/190826/${CHR}_snp --maf 0.${maf} --geno 0.${miss} --make-bed --allow-extra-chr --out ${CHR}_setname_maf${maf}_miss${miss}
  ## 1 of 10k
  gawk 'BEGIN{s=0;chr=""} {if($1==chr){if($4>s+10000){print $2;s=$4}}else{chr=$1;print $2;s=$4}}' ${CHR}_setname_maf${maf}_miss${miss}.bim > ${CHR}_1in10k_maf${maf}_miss${miss}.prune.in
  plink --bfile ${CHR}_setname_maf${maf}_miss${miss} --allow-extra-chr --extract ${CHR}_1in10k_maf${maf}_miss${miss}.prune.in --make-bed --out ${CHR}_1in10k_maf${maf}_miss${miss}_pruned
  plink --bfile ${CHR}_1in10k_maf${maf}_miss${miss}_pruned --distance square 1-ibs flat-missing --out ${CHR}_1in10k_maf${maf}_miss${miss}_pruned --allow-extra-chr
  Rscript ${PLOTSCRIPT} -d ${CHR}_1in10k_maf${maf}_miss${miss}_pruned.mdist -s ${SAMPLEORDER} -a /data2/rawdata2/sample_metadata/withXN/metadata_all_v3_5groups.txt -p ${SAMPLEORDER} -c ${COLORFILE} -P T -o ${CHR}_1in10k_maf${maf}_miss${miss}_pruned -t ${CHR}_1in10k_maf${maf}_miss${miss}_pruned)
done
qpdf --empty --pages {DD,AABBDD,AABB}_1in10k_maf${maf}_miss${miss}_pruned.pdf -- 1in10k_maf${maf}_miss${miss}_pruned.pdf) &
sleep 5
done
done
wait

