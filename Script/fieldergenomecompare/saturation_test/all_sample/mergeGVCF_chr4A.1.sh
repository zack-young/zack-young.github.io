#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T GenotypeGVCFs \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/05_sample/chr4A.1.g.ch.vcf.gz \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/10_sample/chr4A.1.g.ch.vcf.gz \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/25_sample/chr4A.1.g.ch.vcf.gz \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/50_sample/chr4A.1.g.ch.vcf.gz \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/75_sample/chr4A.1.g.ch.vcf.gz \
  --variant /data/user/yangzz/mapping/fieldergenomecompare/saturation_test/85_sample/chr4A.1.g.ch.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/compare/s30/06.callsnp/chr4A.1.g.vcf.gz \
  -o chr4A.1.raw.vcf.gz
