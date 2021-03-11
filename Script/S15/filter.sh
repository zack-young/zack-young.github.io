#!/bin/bash

for vcf in S4/*raw.vcf; do
  ./secondstep ${vcf}
done
#awk '($2 > 541000001) && ($2 < 542000000) && ($5 != "./.") && ($6 != "./.") \&& ($7 != "./.") && ($6 != $5) && ($6 != $7){print $0}' chr2B.raw.bcf.filter
