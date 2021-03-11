#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19

#source /WORK/app/osenv/ln1/set2.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start Callsnp"

WP="/WORK/pp192/"
#
CHR=$1
DIR=$2
#
i=$2
#samtools view -h ${DIR}/${CHR}.dedup.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > ${DIR}/${CHR}.uniq.bam
#samtools index ${DIR}/${CHR}.uniq.bam
samtools flagstat ${DIR}/${CHR}.dedup.bam > ${DIR}/${CHR}.dedup.flagstat
#
java -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups I=${DIR}/${CHR}.dedup.bam O=${DIR}/${CHR}.dedup_ID.bam RGID=${i} RGLB=${i} RGPL=illumina RGSM=${i} RGPU=unit1

samtools index ${DIR}/${CHR}.dedup_ID.bam
echo $CHR
#  -jdk_deflater -jdk_inflater are new features in 3.8 but not support in GZ
java -Xmx40g \
  -Djava.io.tmpdir=/data/user/yangzz/tmp \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -jdk_deflater \
  -jdk_inflater \
  -nct 6 \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T HaplotypeCaller -ERC GVCF -L ${CHR} \
  -I ${DIR}/${CHR}.dedup_ID.bam \
  -o ${DIR}/${CHR}.g.vcf.gz

wait
#
#rm ${DIR}/${CHR}.uniq.bam*
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish Callsnp"

