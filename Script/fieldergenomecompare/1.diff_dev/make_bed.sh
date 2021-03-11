#!/bin/bash

if false;then
bedtools makewindows -g /data/genome/wheat/CS_IWGSC/v1/161010_Chinese_Spring_v1.0_pseudomolecules.fasta.fai -w 1000000 > CS_win1M.bed
gawk '{print > $1".1M.bed"}' CS_win1M.bed

bedtools makewindows -g /data/genome/wheat/CS_IWGSC/v1/161010_Chinese_Spring_v1.0_pseudomolecules.fasta.fai -w 20000000 > CS_win20M.bed
gawk '{print > $1".20M.bed"}' CS_win20M.bed

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  declare -i i=0
  total=$(($(wc -l < ${CHR}.1M.bed)/20))
  for ((i=0;i<=$total;i++));do
    # newname=${CHR}.$(gawk -vi=$i 'BEGIN{print tolower(sprintf("%c",i+65))}').bed
    newname=${CHR}.$(printf "%02d" ${i}).bed
    tail -n +$((i*20+1)) ${CHR}.1M.bed| head -n 20 > ${newname}
  done
done
fi

bedtools makewindows -g /data/genome/wheat/CS_IWGSC/v1/161010_Chinese_Spring_v1.0_pseudomolecules.fasta.fai -w 100000 > CS_win100k.bed
gawk '{print > $1".100k.bed"}' CS_win100k.bed

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  declare -i i=0
  total=$(($(wc -l < ${CHR}.100k.bed)/200))
  for ((i=0;i<=$total;i++));do
    newname=${CHR}.$(printf "%02d" ${i}).100k.bed
    tail -n +$((i*200+1)) ${CHR}.100k.bed| head -n 200 > ${newname}
  done
done
