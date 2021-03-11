#!/usr/bin/env bash
set -euxo pipefail

# for i in {1..7};do
# for j in A B D;do
#   echo chr${i}${j}_setname >> AABBDD.txt
#   echo chr${i}${j}_ld_pruned >> AABBDD_ld_pruned.txt
#   echo chr${i}${j}_gene_pruned >> AABBDD_gene_pruned.txt
#   echo chr${i}${j}_1in10k_pruned >> AABBDD_1in10k_pruned.txt
# done
# done

for i in {1..7};do
for j in A B D;do
for k in _1in10k_pruned _gene_pruned _ld_pruned _setname;do
  if [ $j = "A" ] || [ $j = "B" ];then
    WD=/data2/rawdata2/tree/combineXN_190803/AABB_sub_190817/01.bed
  else
    WD=/data2/rawdata2/tree/combineXN_190803/D_sub_190809/01.bed
  fi

  plink --bfile ${WD}/chr${i}${j}${k} --keep AABBDD_keepfile.txt --make-bed --allow-extra-chr --out chr${i}${j}${k}
done
done
done

for i in AABBDD;do
  plink --merge-list ${i}.txt --make-bed --out ${i}_snp --allow-extra-chr &
  plink --merge-list ${i}_ld_pruned.txt --make-bed --out ${i}_snp_ld_pruned --allow-extra-chr &
  plink --merge-list ${i}_gene_pruned.txt --make-bed --out ${i}_snp_gene_pruned --allow-extra-chr &
  plink --merge-list ${i}_1in10k_pruned.txt --make-bed --out ${i}_snp_1in10k_pruned --allow-extra-chr &
done
#
wait
