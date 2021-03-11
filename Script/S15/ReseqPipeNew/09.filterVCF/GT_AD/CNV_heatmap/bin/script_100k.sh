#!/usr/bin/env usr/sh
set -euxo pipefail

# sum bin to 100k
SM=$1
#mkdir $SM
for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  gawk '{sum+=$4} (NR%100)==0{print sum/100; sum=0} END{n=(NR%100);if(n!=0){print sum/n}}' ../readDepth_ori_1k/$SM/${chr}.1.1k.DP > $SM/${chr}.1.100k.tmp
  gawk '{sum+=$4} (NR%100)==0{print sum/100; sum=0} END{n=(NR%100);if(n!=0){print sum/n}}' ../readDepth_ori_1k/$SM/${chr}.2.1k.DP > $SM/${chr}.2.100k.tmp
done

# normlize DP coverage by sample
sum=`gawk '{sum+=$1} END{print sum}' ${SM}/*100k.tmp`
line=`wc -l ${SM}/*100k.tmp|tail -n 1|gawk '{print $1}'`
ave=$(perl -e "print $sum/$line")

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${SM}/${CHR}.1.100k.tmp > ${SM}/${CHR}.100k.norm
  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${SM}/${CHR}.2.100k.tmp >> ${SM}/${CHR}.100k.norm
  sed -i '1i '$SM'' ${SM}/${CHR}.100k.norm
  sed '1 s/^.*$/'$CHR'/' ${SM}/${CHR}.100k.norm > ${SM}/${CHR}.100k.norm.chr
done

## join DP by chr
#while read line;
#  do paste ../*/${line}.100k.norm > ${line}.100k.join; 
#done < chrlst
#mv ${line}.100k.join picture_100k

## join DP info by sample for fingerprint 
paste ${SM}/*.100k.norm.chr > ${SM}/${SM}.100k.norm.chr; cp ${SM}/${SM}.100k.norm.chr picture_100k/${SM}.chr
