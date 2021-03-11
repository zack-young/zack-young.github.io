#!/usr/bin/env usr/sh
set -euxo pipefail

# sum bin to 100k
SM=$1
path_SM="/data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/${1}_${2}"
#mkdir $path_SM
for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
#  gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' /data2/rawdata2/readDepth_ori_1k/$SM/${chr}.1.1k.DP > $SM/${chr}.1.1M.tmp
#  gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' /data2/rawdata2/readDepth_ori_1k/$SM/${chr}.2.1k.DP > $SM/${chr}.2.1M.tmp
    gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' $path_SM/${chr}.1.1k.DP > $path_SM/${chr}.1.1M.tmp
    gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' $path_SM/${chr}.2.1k.DP > $path_SM/${chr}.2.1M.tmp
done

# normlize DP coverage by sample
sum=`gawk '{sum+=$1} END{print sum}' ${path_SM}/*1M.tmp`
line=`wc -l ${path_SM}/*1M.tmp|tail -n 1|gawk '{print $1}'`
ave=$(perl -e "print $sum/$line")

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${path_SM}/${CHR}.1.1M.tmp > ${path_SM}/${CHR}.1M.norm
  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${path_SM}/${CHR}.2.1M.tmp >> ${path_SM}/${CHR}.1M.norm
  sed -i '1i '$SM'' ${path_SM}/${CHR}.1M.norm
  sed '1 s/^.*$/'$CHR'/' ${path_SM}/${CHR}.1M.norm > ${path_SM}/${CHR}.1M.norm.chr
done

## join DP by chr
#while read line;
#  do paste ../*/${line}.100k.norm > ${line}.100k.join; 
#done < chrlst
#mv ${line}.100k.join picture_100k

paste ${path_SM}/*.1M.norm.chr > ${path_SM}/${SM}.1M.norm.chr
#cp ${SM}/${SM}.1M.norm.chr picture_100k/${SM}.chr
