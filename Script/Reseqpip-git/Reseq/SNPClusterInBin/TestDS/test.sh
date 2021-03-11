#!/usr/bin/env sh


# test SNPClusterInBin.py
../SNPClusterInBin.py -b 1000 -C 0.001 -i test.vcf > tmp_output.SNPCluster
cat test.vcf | ../SNPClusterInBin.py -b 1000 -C 0.05 -o tmp_output.SNPCluster.gz


rm tmp_*

