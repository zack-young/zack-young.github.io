#!/bin/bash

if true;then
#for CHR in chr1A chr2A chr3A chr4A chr5A chr6A chr7A;do
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#for CHR in chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
  #n=$(wc -l < ${CHR}.1M.bed)
  #parallel -j procfile bash ./script.sh ${CHR} {} ::: $(eval echo {1..$n})
  python ./concat_dist.py -c ${CHR} -b 1 &
  sleep 10
done
fi

