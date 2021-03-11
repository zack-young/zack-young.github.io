#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T GenotypeGVCFs \
  --variant /data2/public/ReseqData/workflowFile/YM06/06.callsnp/chr7A.2.g.vcf.gz \
  --variant /data2/public/ReseqData/workflowFile/YM07/06.callsnp/chr7A.2.g.vcf.gz \
  -o chr7A.2.raw.vcf.gz
