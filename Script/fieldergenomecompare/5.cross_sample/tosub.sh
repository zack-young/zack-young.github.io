#!/bin/bash

for CHR in chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
  python ./concat_dist.py -c ${CHR} &
  sleep 10
done
wait

# convert dist between bins to dist between samples.
# python ./concat_dist_to_SM.py
