#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T GenotypeGVCFs \
  --variant /data/user/yangzz/mapping/Reseq_data/CP12/06.callsnp/chr5A.2.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SN224/06.callsnp/chr5A.2.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/NZ10/06.callsnp/chr5A.2.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S136/06.callsnp/chr5A.2.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/TE20/06.callsnp/chr5A.2.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SNFu63/06.callsnp/chr5A.2.g.vcf.gz \
  -o chr5A.2.raw.vcf.gz
