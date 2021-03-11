#!/bin/bash
set -euo pipefail
WD="/data2/rawdata2/bin/datamash-1.4"
if true;then
#gawk '{a[$1]} END{for(i in a){for(j in a){if(i!=j){print i"\t"j}}}}' D_WGS_GD_linen.txt|shuf -n 1000 > D_WGS_GD_linen_rand1000.txt
#gawk 'ARGIND==1{a[$1]} ARGIND==2{b[$1]} END{for(i in a){for(j in b){if(i!=j){print i"\t"j}}}}' D_WGS_GH_linen.txt D_WGS_GD_linen.txt |shuf -n 1000 > D_WGS_GHD_linen_rand1000.txt

for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1D chr2D chr3D chr4D chr5D chr6D chr7D chr1B chr2B chr3B chr4B chr5B chr6B chr7B;do
  (for i in ../raw/${CHR}*1M.dist.txt;do
    gawk -vchr=$CHR 'ARGIND==1{a[$1"."$2]=NR} ARGIND==2{for(i=1;i<=NF;i++){if(FNR"."i in a){print int($i)"\t"a[FNR"."i] >> chr".dist_withsm.txt"}}}' 10_near_linen.txt $i
done ) &
done

wait

cat chr?A.dist_withsm.txt|gawk '{print $1 > "A_pair"$2}'
cat chr?B.dist_withsm.txt|gawk '{print $1 > "B_pair"$2}'
cat chr?D.dist_withsm.txt|gawk '{print $1 > "D_pair"$2}'
fi

if false;then
for i in {1..1247};do
  cat A_pair${i} | ${WD}/datamash --sort --group 1 count 1 |csvtk join -H -f "1;1" -t -k --na 0 group.txt -| cut -f 2 | ${WD}/datamash transpose
done > A_sample_by_count.txt


for i in {1..1247};do
  cat B_pair${i} | ${WD}/datamash --sort --group 1 count 1 |csvtk join -H -f "1;1" -t -k --na 0 group.txt -| cut -f 2 | ${WD}/datamash transpose
done > B_sample_by_count.txt

for i in {1..1247};do
  cat D_pair${i} | ${WD}/datamash --sort --group 1 count 1 |csvtk join -H -f "1;1" -t -k --na 0 group.txt -| cut -f 2 | ${WD}/datamash transpose
done > D_sample_by_count.txt
fi
