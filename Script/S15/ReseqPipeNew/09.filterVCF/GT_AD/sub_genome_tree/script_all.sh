#!/usr/bin/env bash
set -euxo pipefail

PLOTSCRIPT=/data2/rawdata2/tree/script/NJtree_iTOL.R
COLORFILE=/data2/rawdata2/tree/combineXN_190803/colorF_5cols.txt
SAMPLEORDER=/data2/rawdata2/tree/combineXN_190803/D_sub_190809/03.plot/sample_order_hex.txt

for maf in 003;do
for miss in 05 1 2 3 4 5;do
(for CHR in AABB DD AABBDD;do
  (plink --bfile /data2/rawdata2/tree/combineXN_190803/rawdata/190826/${CHR}_snp --maf 0.${maf} --geno 0.${miss} --make-bed --allow-extra-chr --out ${CHR}_setname_maf${maf}_miss${miss}
  plink --bfile ${CHR}_setname_maf${maf}_miss${miss} --distance square 1-ibs flat-missing --out ${CHR}_setname_maf${maf}_miss${miss} --allow-extra-chr
  Rscript ${PLOTSCRIPT} -d ${CHR}_setname_maf${maf}_miss${miss}.mdist -s ${SAMPLEORDER} -a /data2/rawdata2/sample_metadata/withXN/metadata_all_v3_5groups.txt -p ${SAMPLEORDER} -c ${COLORFILE} -P T -o ${CHR}_setname_maf${maf}_miss${miss} -t ${CHR}_setname_maf${maf}_miss${miss} )
done
qpdf --empty --pages {DD,AABBDD,AABB}_setname_maf${maf}_miss${miss}.pdf -- all_maf${maf}_miss${miss}.pdf) &
wait_all
done
done
wait

