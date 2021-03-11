#!/usr/bin/env bash
set -euxo pipefail

WD=/data2/rawdata2/tree/combineXN_190803/AABBDD_sub_190817/

for i in AABBDD;do
  for j in snp_1in10k_pruned snp snp_gene_pruned snp_ld_pruned;do
    plink --bfile ${WD}/01.bed/${i}_${j} --distance square 1-ibs flat-missing --allow-extra-chr --out ${i}_${j} &
  done
done

wait
