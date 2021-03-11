#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# Guo, Weilong; guoweilong@126.com; 2017-10-19

#source /WORK/app/osenv/ln1/set2.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start Callsnp"

WP="/data/user/yangzz/mapping/Reseqpip-git/Reseq/ReseqPipeNew"
#
CHR=$1
set +e
mkdir $2
set -e
cd $2
#
samtools view -h ${WP}/05.*/${2}/${CHR}.dedup.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > ${CHR}.uniq.bam
samtools index ${CHR}.uniq.bam
samtools flagstat ${CHR}.uniq.bam > ${CHR}.flagstat
#
echo $CHR
#  -jdk_deflater -jdk_inflater are new features in 3.8 but not support in GZ
java -Xmx60g \
  -Djava.io.tmpdir=${WP}/05.*/${2}/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -jdk_deflater \
  -jdk_inflater \
  -nct 24 \
  -R ${WP}/06.*/iwgsc_refseqv2.0_all_chromosomes_parts_sed.fa \
  -T HaplotypeCaller -ERC GVCF -L ${CHR} \
  -I ${WP}/05.*/$2/${CHR}.dedup.bam \
  -o ${CHR}.g.vcf.gz
#
wait
#
rm ${CHR}.uniq.bam*
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish Callsnp"

