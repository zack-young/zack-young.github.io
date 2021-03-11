#!/usr/bin/env bash
#
WP=/data/user/yangzz/mapping/Reseq_data
#/data2/rawdata2/workflowFile
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start merge gVCF"
#
batch=$1
SMlst=$2
START=$3
END=$4
#
for CHR in `echo ${batch} | tr "," " "`; do
    ((echo '#!/usr/bin/env sh';
    echo "java  \\";
    echo "  -Djava.io.tmpdir=/home/wangzh/tmp/ \\";
    echo "  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \\";
    echo "  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \\";
    echo "  -L ${CHR}:${START}-${END} \\";
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
#echo "  --variant /data2/public/ReseqData/workflowFile/Ak58Simu/06.callsnp/${CHR}.g.vcf.gz \\" >> mergeGVCF_${CHR}.sh
done
#
wait
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish merge gVCF"

# how to use variable as variable name ?
# eval batch='$'$i
