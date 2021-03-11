#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T GenotypeGVCFs \
  --variant /data/user/yangzz/mapping/Reseq_data/NZ01/06.callsnp/chr4B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/XM01/06.callsnp/chr4B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766495/06.callsnp/chr4B.1.g.vcf.gz \
  -o chr4B.1.raw.vcf.gz
