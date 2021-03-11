#!/usr/bin/env bash
#set -euxo pipefail

line=$1
CHR=$2
mkdir /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/$line

# calculate DP coverage for each sample & chr
bedtools coverage -a ~/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/bed/${CHR}.1k.bed -b ~/mapping/Reseqpip-git/Reseq/ReseqPipeNew/05.mergeAsplit/${line}/${CHR}.*.bam -counts > /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/${line}/${CHR}.1k.DP
