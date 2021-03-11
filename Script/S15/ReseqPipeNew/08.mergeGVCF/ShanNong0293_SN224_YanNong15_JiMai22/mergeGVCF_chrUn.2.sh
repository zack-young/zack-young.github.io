#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T GenotypeGVCFs \
  --variant /data2/public/ReseqData/workflowFile/CP12/06.callsnp/chrUn.2.g.vcf.gz \
  --variant /data2/public/ReseqData/workflowFile/SN224/06.callsnp/chrUn.2.g.vcf.gz \
  --variant /data2/public/ReseqData/workflowFile/NZ10/06.callsnp/chrUn.2.g.vcf.gz \
  --variant /data2/public/ReseqData/workflowFile/S136/06.callsnp/chrUn.2.g.vcf.gz \
  -o chrUn.2.raw.vcf.gz
