#!/usr/bin/env usr/sh
set -euxo pipefail

CHR=$1; echo "CHR : ${CHR}"
SM=${CHR}.raw.vcf.gz; echo "Sample: ${SM}"

# SNP
java -Djava.io.tmpdir=/home/wangzh/tmp/ \
     -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
     -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
     -T SelectVariants \
     -V ${SM} \
     -selectType SNP \
     -o ${CHR}.snp.vcf.gz

java -Djava.io.tmpdir=/home/wangzh/tmp/ \
     -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
     -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
     -T VariantFiltration \
     -V ${CHR}.snp.vcf.gz \
     --logging_level ERROR \
     --filterExpression " QD < 2.0 " \
     --filterName "QD_filter" \
     --filterExpression " FS > 60.0 " \
     --filterName "FS_filter" \
     --filterExpression " MQRankSum < -12.5 " \
     --filterName "MQRankSum_filter" \
     --filterExpression " ReadPosRankSum < -8.0 " \
     --filterName "ReadPosRankSum_filter" \
     --filterExpression " SOR > 3.0 " \
     --filterName "SOR_filter" \
     --filterExpression " MQ < 40.0 " \
     --filterName "MQ_filter" \
     --genotypeFilterExpression " DP > 30 || DP < 3 " \
     --genotypeFilterName "DP_filter" \
     --setFilteredGtToNocall \
     --clusterSize 3 \
     --clusterWindowSize 10 \
     -o ${CHR}.snp.filter.vcf.gz

java -Djava.io.tmpdir=/home/wangzh/tmp/ \
     -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
     -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
     -T SelectVariants \
     -V ${CHR}.snp.filter.vcf.gz \
     --excludeFiltered \
     --excludeNonVariants \
     --removeUnusedAlternates \
     -o ${CHR}.snp.filter.final.vcf.gz

rm ${CHR}.snp.vcf* ${CHR}.snp.filter.vcf*

# INDEL
java -Djava.io.tmpdir=/home/wangzh/tmp/ \
     -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
     -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
     -T SelectVariants \
     -V ${SM} \
     -selectType INDEL \
     -o ${CHR}.indel.vcf.gz

java -Djava.io.tmpdir=/home/wangzh/tmp/ \
     -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
     -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
     -T VariantFiltration \
     -V ${CHR}.indel.vcf.gz \
     --logging_level ERROR \
     --filterExpression " QD < 2.0  " \
     --filterName "QD_filter" \
     --filterExpression " FS > 200.0 " \
     --filterName "FS_filter" \
     --filterExpression " ReadPosRankSum < -20.0 " \
     --filterName "ReadPosRankSum_filter" \
     --genotypeFilterExpression " DP > 30 || DP < 3 " \
     --genotypeFilterName "DP_filter" \
     --setFilteredGtToNocall \
     -o ${CHR}.indel.filter.vcf.gz

java -Djava.io.tmpdir=/home/wangzh/tmp/ \
     -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
     -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
     -T SelectVariants \
     -V ${CHR}.indel.filter.vcf.gz \
     --excludeFiltered \
     --excludeNonVariants \
     --removeUnusedAlternates \
     -o ${CHR}.indel.filter.final.vcf.gz

rm ${CHR}.indel.vcf* ${CHR}.indel.filter.vcf*