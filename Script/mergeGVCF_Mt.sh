#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T GenotypeGVCFs \
  --variant /data/rawdata/workflowFile/LX987/06.callsnp/Mt.g.vcf.gz \
  --variant /data/rawdata/workflowFile/s3097-1/06.callsnp/Mt.g.vcf.gz \
  --variant /data/rawdata/workflowFile/S140/06.callsnp/Mt.g.vcf.gz \
  -o Mt.raw.vcf.gz
