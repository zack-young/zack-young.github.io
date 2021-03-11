#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T GenotypeGVCFs \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/05_sample/chr7A.2.g.ch.vcf.gz \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/10_sample/chr7A.2.g.ch.vcf.gz \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/25_sample/chr7A.2.g.ch.vcf.gz \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/50_sample/chr7A.2.g.ch.vcf.gz \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/75_sample/chr7A.2.g.ch.vcf.gz \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/85_sample/chr7A.2.g.ch.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/compare/s30/06.callsnp/chr7A.2.g.vcf.gz \
  -o chr7A.2.raw.vcf.gz
