#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

WP="/WORK/pp192/"
source /WORK/app/osenv/ln1/set2.sh

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start SplitGVCF"

for VCF in ../06.callsnp/*.vcf.gz;do
    (CHR=`basename $VCF|sed s/.g.vcf.gz//g`
    python SplitGvcfByBin.py -i ../06.callsnp/${CHR}*.vcf.gz -B 20000000;
    for ID in ${CHR}.??.g.vcf;do
        bgzip ${ID};
        tabix -p vcf ${ID}.gz
    done ) &
done
#
wait

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Finish SplitGVCF"