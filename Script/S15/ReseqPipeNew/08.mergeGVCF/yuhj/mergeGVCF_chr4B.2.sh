#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T GenotypeGVCFs \
  --variant /data/user/yangzz/mapping/Reseq_data/PH154/06.callsnp/chr4B.2.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH06/06.callsnp/chr4B.2.g.vcf.gz \
  -o chr4B.2.raw.vcf.gz
