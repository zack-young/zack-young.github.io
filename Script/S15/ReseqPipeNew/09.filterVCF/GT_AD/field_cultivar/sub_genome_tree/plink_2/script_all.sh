#!/usr/bin/env bash
set -euxo pipefail

PLOTSCRIPT=/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/sub_genome_tree/NJtree_iTOL.R
COLORFILE=/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/sub_genome_tree/coloryzz_6cols.txt
SAMPLEORDER=/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/cultivar_list

for maf in 003;do
for miss in 05 1 2 3 4 5;do
(for CHR in AABB DD AABBDD;do
  (/data/user/yangzz/worktools/Plink/plink --bfile /data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/sub_genome_tree/plink_1/${CHR}_snp --maf 0.${maf} --geno 0.${miss} --make-bed --allow-extra-chr --out ${CHR}_setname_maf${maf}_miss${miss}
  /data/user/yangzz/worktools/Plink/plink --bfile ${CHR}_setname_maf${maf}_miss${miss} --distance square 1-ibs flat-missing --out ${CHR}_setname_maf${maf}_miss${miss} --allow-extra-chr
  Rscript ${PLOTSCRIPT} -d ${CHR}_setname_maf${maf}_miss${miss}.mdist -s ${SAMPLEORDER} -a /data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/metadata_cultivar_nogroup.txt -p ${SAMPLEORDER} -c ${COLORFILE} -P T -o ${CHR}_setname_maf${maf}_miss${miss} -t ${CHR}_setname_maf${maf}_miss${miss} )
done
qpdf --empty --pages {DD,AABBDD,AABB}_setname_maf${maf}_miss${miss}.pdf -- all_maf${maf}_miss${miss}.pdf) &
wait_all
done
done
wait
