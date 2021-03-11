#!/usr/bin/env bash
set -euxo pipefail

PLOTSCRIPT=/data2/rawdata2/tree/script/NJtree_iTOL.R
COLORFILE=/data2/rawdata2/tree/combineXN_190803/colorF_5cols.txt
SAMPLEORDER=/data2/rawdata2/tree/combineXN_190803/D_sub_190809/03.plot/sample_order_hex.txt

for maf in 01 02 03 04 05 1;do
for miss in 05 1 2 3 4 5;do
(for CHR in AABB DD AABBDD;do
  (#plink --bfile /data2/rawdata2/tree/combineXN_190803/rawdata/190826/${CHR}_snp --maf 0.${maf} --geno 0.${miss} --make-bed --allow-extra-chr --out ${CHR}_setname_maf${maf}_miss${miss}
  ## gene
  bedtools intersect -b /data2/rawdata2/database/CSv1p1.all.genes.merged.bed -a <(gawk -vOFS="\t" '{print $1,$4,$4}' ${CHR}_setname_maf${maf}_miss${miss}.bim) -sorted | gawk '{print $1":"$2}' > ${CHR}_gene_maf${maf}_miss${miss}.prune.in
  plink --bfile ${CHR}_setname_maf${maf}_miss${miss} --allow-extra-chr --extract ${CHR}_gene_maf${maf}_miss${miss}.prune.in --make-bed --out ${CHR}_gene_maf${maf}_miss${miss}_pruned
  plink --bfile ${CHR}_gene_maf${maf}_miss${miss}_pruned --distance square 1-ibs flat-missing --out ${CHR}_gene_maf${maf}_miss${miss}_pruned --allow-extra-chr
  Rscript ${PLOTSCRIPT} -d ${CHR}_gene_maf${maf}_miss${miss}_pruned.mdist -s ${SAMPLEORDER} -a /data2/rawdata2/sample_metadata/withXN/metadata_all_v3_5groups.txt -p ${SAMPLEORDER} -c ${COLORFILE} -P T -o ${CHR}_gene_maf${maf}_miss${miss}_pruned -t ${CHR}_gene_maf${maf}_miss${miss}_pruned)
done
qpdf --empty --pages {DD,AABBDD,AABB}_gene_maf${maf}_miss${miss}_pruned.pdf -- gene_maf${maf}_miss${miss}_pruned.pdf) &
sleep 5
done
done
wait

