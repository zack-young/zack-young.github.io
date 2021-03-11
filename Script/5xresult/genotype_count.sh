#!/usr/bin/bash
for item in *.raw.bcf.filter; do
    CHR=`basename $item | sed s/.raw.bcf.filter//g`
    echo ${CHR}'_threediff' >> genotype_count.txt
    less ${item}|gawk -vOFS="\t" '$5 != $6 && $6 != $7 && $5 != $7'|gawk -vOFS="\t" '$5 != "./."&&$6 != "./."&&$7 != "./."{print}' | wc -l >> genotype_count.txt
    echo ${CHR}'_twodiff' >> genotype_count.txt
    less ${item}|gawk -vOFS="\t" '$5 != $6 && $6 != $7'|gawk -vOFS="\t" '$5 != "./."&&$6 != "./."&&$7 != "./."{print}' | wc -l >> genotype_count.txt
    echo ${CHR}'_line' >> genotype_count.txt
    less ${item}| wc -l >> genotype_count.txt
done
