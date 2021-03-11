#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T GenotypeGVCFs \
  --variant /data2/public/ReseqData/workflowFile/NZ33/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/public/ReseqData/workflowFile/S64/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/public/ReseqData/workflowFile/NZ34/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/public/ReseqData/workflowFile/NZ18/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/public/ReseqData/workflowFile/NZ19/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/Nongke145/SRR10766582/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/Nongke145/SRR10766583/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/Nongke145/SRR10766588/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/Nongke145/SRR10766624/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/Nongke145/SRR10766512/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/Nongke145/SRR10766529/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/Nongke145/SRR10766555/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/Nongke145/SRR10766564/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/Nongke145/SRR10766628/06.callsnp/chr2A.2.g.vcf.gz \
  --variant /data2/user2/xiexm/pro_ject/resequencing/Nongke145/SRR10766631/06.callsnp/chr2A.2.g.vcf.gz \
  -o chr2A.2.raw.vcf.gz
