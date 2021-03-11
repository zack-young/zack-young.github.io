#!/bin/bash

PLOTSCRIPT=/data2/rawdata2/tree/script/NJtree_iTOL.R
COLORFILE=/data2/rawdata2/tree/combineXN_190803/colorF_4cols.txt
CHR=AABBDD

Rscript ${PLOTSCRIPT} -d ../02.dist/${CHR}_snp_gene_pruned.mdist -s sample_order.txt -a /data2/rawdata2/sample_metadata/withXN/metadata_all_v3.txt -p sample_order.txt -c ${COLORFILE} -P T -o ${CHR}_snp_gene_pruned

Rscript ${PLOTSCRIPT} -d ../02.dist/${CHR}_snp_1in10k_pruned.mdist -s sample_order.txt -a /data2/rawdata2/sample_metadata/withXN/metadata_all_v3.txt -p sample_order.txt -c ${COLORFILE} -P T -o ${CHR}_snp_1in10k_pruned

Rscript ${PLOTSCRIPT} -d ../02.dist/${CHR}_snp_ld_pruned.mdist -s sample_order.txt -a /data2/rawdata2/sample_metadata/withXN/metadata_all_v3.txt -p sample_order.txt -c ${COLORFILE} -P T -o ${CHR}_snp_ld_pruned

Rscript ${PLOTSCRIPT} -d ../02.dist/${CHR}_snp.mdist -s sample_order.txt -a /data2/rawdata2/sample_metadata/withXN/metadata_all_v3.txt -p sample_order.txt -c ${COLORFILE} -P T -o ${CHR}_snp
