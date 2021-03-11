#!/usr/bin/env bash
set -euxo pipefail

CHR=$1
FILE=${CHR}.bcf.gz

#echo ${2:-IWGSC}

bcftools view $FILE --threads 4 \
    | java -jar  /home/wangzh/bin/snpEff_latest_core/snpEff/snpEff.jar -no-downstream -no-upstream  ${2:-IWGSC} - -noStats \
    | bcftools view --threads 4 -o ${CHR}.ann.bcf.gz -Ob
#| java -jar  /home/wangzh/bin/snpEff_latest_core/snpEff/snpEff.jar -no-downstream -no-upstream -t ${2:-IWGSC} - -csvStats ${CHR}.ann.vcf.csv  \
