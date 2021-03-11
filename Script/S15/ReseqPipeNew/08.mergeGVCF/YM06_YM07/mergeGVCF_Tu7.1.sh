#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /data2/rawdata2/genome/Tu/Tu_split.fa \
  -T GenotypeGVCFs \
  --variant /data2/public/ReseqData/workflowFile/PH11_to_Tu/06.callsnp/Tu7.1.g.vcf.gz \
  --variant /data2/public/ReseqData/workflowFile/PH12_to_Tu/06.callsnp/Tu7.1.g.vcf.gz \
  --variant /data2/public/ReseqData/workflowFile/PH15_to_Tu/06.callsnp/Tu7.1.g.vcf.gz \
  -o Tu7.1.raw.vcf.gz
