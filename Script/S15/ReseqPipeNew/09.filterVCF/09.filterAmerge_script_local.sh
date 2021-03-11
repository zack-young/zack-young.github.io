#!/usr/bin/env bash
#set -euxo pipefail


echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start filter and merge GVCF"

#    echo '#SBATCH -N 1' >> filterGVCF_${CHR}.sh
#    echo '#SBATCH -n 1' >> filterGVCF_${CHR}.sh
#    echo '#SBATCH -c 24' >> filterGVCF_${CHR}.sh
#    echo '#SBATCH -t 2400' >> filterGVCF_${CHR}.sh
#    echo 'WP="/WORK/pp192/"' >> filterGVCF_${CHR}.sh
#for GVCF in `ls ../08.mergeGVCF/*.raw.vcf.gz`;do
for GVCF in chr3A.1.raw.vcf.gz chr3B.1.raw.vcf.gz chr4B.2.raw.vcf.gz chr4D.1.raw.vcf.gz chr7A.2.raw.vcf.gz;do
    CHR=`basename $GVCF|sed s/.raw.vcf.gz//g`
    echo ${CHR}
    sh VCFfilter_local.sh ${CHR}> ${CHR}.o 2>${CHR}.e &
    wait_all
#     yhbatch -n 1 -e ${CHR}.e -o ${CHR}.o ./filterGVCF_${CHR}.sh
        #`
        #sed "s/\$CHR/$CHR/g;s/\$ID/$ID/g" VCFfilter.sh > VCFfilter_${CHR}_${ID}.sh
        #cat VCFfilter_${CHR}_${ID}.sh >> filterGVCF_${CHR}.sh
        #rm VCFfilter_${CHR}_${ID}.sh
done
    #echo "wait"  >> filterGVCF_${CHR}.sh
wait
#
# wait for all done
# yhbatch -n 1 concatGVCF.sh -e concat.e -o concat.o ./concatGVCF.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] End filter and merge GVCF"
