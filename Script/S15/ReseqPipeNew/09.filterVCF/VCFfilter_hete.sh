#!/usr/bin/env sh
set -euxo pipefail

#CHR=$1; echo "CHR : ${CHR}"
#SM=${CHR}.raw.vcf.gz; echo "Sample: ${SM}"

# SNP
java -Djava.io.tmpdir=/home/wangzh/tmp/ \
     -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
     -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
     -T SelectVariants \
     -V chr1A.1.snp.filter.final.vcf.gz \
     -select 'vc.getGenotype("S1").getAD().0 == 1 || vc.getGenotype("S1").getAD().1 == 1' \
     -o chr1A.1.snp.hete.vcf.gz

