#!/usr/bin/env bash

CHRl=$1
for CHR in `echo $CHRl|tr "," " "`;do
    yhbatch -n 1 -e concat_${CHR}.e -o concat_${CHR}.o ./concatVCF.sh ${CHR}
done