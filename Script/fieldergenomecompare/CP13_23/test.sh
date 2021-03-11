#!/usr/bin/env usr/sh
set -euxo pipefail

# sum bin to 100k
SM=$1
for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
    gawk '{sum+=$4} (NR%10)==0{print sum/10; sum=0} END{n=(NR%10);if(n!=0){print sum/n}}' ${SM}_2020-04-17162314/${chr}.1.1M.tmp > ${SM}_2020-04-17162314/${chr}.1.10M.tmp
    gawk '{sum+=$4} (NR%10)==0{print sum/10; sum=0} END{n=(NR%10);if(n!=0){print sum/n}}' ${SM}_2020-04-17162314/${chr}.2.1M.tmp > ${SM}_2020-04-17162314/${chr}.2.10M.tmp
#    gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' $path_SM/${chr}.1.1k.DP > ${1}/${chr}.1.1M.tmp
#    gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' $path_SM/${chr}.2.1k.DP > ${1}/${chr}.2.1M.tmp
    #cat ${1}/${chr}.1.1M.tmp ${1}/${chr}.2.1M.tmp > ${1}/${chr}.1M.DP
done
cat ${SM}_2020-04-17162314/chr*.10M.tmp > ${SM}_2020-04-17162314/combine_10M_DP
ave=`./count_CNV.py --normalize on -i ${SM}_2020-04-17162314/combine_10M_DP`
# normlize DP coverage by sample
#sum=`gawk '{sum+=$1} END{print sum}' ${path_SM}/*1k.DP`
#line=`awk '{if ($4!=0) print $0}' ${path_SM}/chr*.1k.DP | wc -l |tail -n 1|gawk '{print $1}'`
#ave=$(perl -e "print $sum/$line")

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
#for CHR in chr1B;do
  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${SM}_2020-04-17162314/${CHR}.1.10M.tmp > ${SM}_2020-04-17162314/${CHR}.10M.norm
  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${SM}_2020-04-17162314/${CHR}.2.10M.tmp >> ${SM}_2020-04-17162314/${CHR}.10M.norm
  sed -i '1i '$SM'' ${SM}_2020-04-17162314/${CHR}.10M.norm
  sed '1 s/^.*$/'$CHR'/' ${SM}_2020-04-17162314/${CHR}.10M.norm > ${SM}_2020-04-17162314/${CHR}.10M.norm.chr
done

## join DP by chr
#while read line;
#  do paste ../*/${line}.100k.norm > ${line}.100k.join; 
#done < chrlst
#mv ${line}.100k.join picture_100k

#paste ${path_SM}/*.1M.norm.chr > ${path_SM}/${SM}.1M.norm.chr
#cp ${SM}/${SM}.1M.norm.chr picture_100k/${SM}.chr
