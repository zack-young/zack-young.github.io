#!/usr/bin/env bash
set -euxo pipefail

CHR=$1
FILE=${CHR}.bcf.gz

bcftools view $FILE --threads 4 \
    | java -jar /home/wangzh/bin/snpEff_latest_core/snpEff/snpEff.jar -t ${2:-IWGSC} - \
    | bcftools view --threads 4 -o ${CHR}.ann.bcf.gz -Ob