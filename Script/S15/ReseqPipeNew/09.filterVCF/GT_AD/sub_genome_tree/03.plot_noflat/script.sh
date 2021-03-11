#!/bin/bash
set -x

PLOTSCRIPT=/data2/rawdata2/tree/script/NJtree_iTOL.R
COLORFILE=/data2/rawdata2/tree/combineXN_190803/colorF_4cols.txt

for CHR in AABBDD;do
  for j in snp_1in10k_pruned snp snp_gene_pruned snp_ld_pruned;do
    Rscript ${PLOTSCRIPT} -d ../02.dist_noflat/${CHR}_${j}.mdist -s sample_order.txt -a /data2/rawdata2/sample_metadata/withXN/metadata_all_v3.txt -p sample_order.txt -c ${COLORFILE}  -P T -o ${CHR}_${j}
  done
done
