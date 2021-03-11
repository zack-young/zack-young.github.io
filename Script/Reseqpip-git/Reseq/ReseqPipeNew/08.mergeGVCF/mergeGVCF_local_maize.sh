#!/usr/bin/env bash
#
WP=/data/rawdata/workflowFile/noneWheat/
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start merge gVCF"
#
batch=$1
SMlst=$2
#
for CHR in `echo ${batch} | tr "," " "`; do
    ((echo '#!/usr/bin/env sh';
    echo "java  \\";
    echo "  -Djava.io.tmpdir=/home/wangzh/tmp/ \\";
    echo "  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \\";
    echo "  -R /data/genome/maize/B73/AGPv4/Zea_mays.AGPv4.dna.toplevel.fa \\";
    echo "  -T GenotypeGVCFs \\";) > mergeGVCF_${CHR}.sh
    #
    (for SM in `echo $SMlst | tr "," " "`; do
    echo "  --variant ${WP}/${SM}/06.callsnp/${CHR}.g.vcf.gz \\";
    done) >> mergeGVCF_${CHR}.sh
    #
    echo "  -o ${CHR}.raw.vcf.gz"  >> mergeGVCF_${CHR}.sh
    #
    sh mergeGVCF_${CHR}.sh > ${CHR}.o 2>${CHR}.e) &
    wait_all
done
#
wait
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish merge gVCF"