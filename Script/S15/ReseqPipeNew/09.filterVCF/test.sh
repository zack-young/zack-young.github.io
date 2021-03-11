#!/usr/bin/env sh

#CHR=$1; echo "CHR : ${CHR}"
#SM=$2.raw.vcf.gz; echo "Sample: ${SM}"
#for GVCF in chr3A.1.raw.vcf.gz chr3B.1.raw.vcf.gz chr4B.2.raw.vcf.gz chr4D.1.raw.vcf.gz chr7A.2.raw.vcf.gz;do
#    CHR=`basename $GVCF|sed s/.raw.vcf.gz//g`
#    echo ${CHR}
#done

mkdir raw$1.txt
