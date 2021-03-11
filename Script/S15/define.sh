#!/bin/bash
#./mofaso_third.py -b 1000000 -i /data/user/yangzz/mapping/s1_4/S1/chr1A.raw.vcf.filter  -s 5 -o tes11
for VCF in /data/user/yangzz/mapping/s1_4/S3/chr??.raw.vcf.filter; do
  CHR=`echo $VCF | sed s/.raw.vcf.filter//g`
  echo "${VCF} ${CHR}"
  ./mofaso_hete.py -b 1000000 -i ${VCF} -s 5 -o ${CHR}_count
done

