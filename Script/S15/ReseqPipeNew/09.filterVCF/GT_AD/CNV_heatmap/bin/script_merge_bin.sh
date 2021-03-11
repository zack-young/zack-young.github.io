#!/bin/bash
set -euxo pipefail

bin=5k
bin1=$(echo ${bin} | sed 's/g/ * 1000 m/;s/m/ * 1000 k/;s/k/ * 1000/;'|bc)
if true;then
for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  # to formated bed
  #paste ../CS_to_CS/${CHR}.${bin}.norm ../Zang1817_to_CS/${CHR}.${bin}.norm | gawk -vOFS="\t" -vc=$CHR -vb=$bin1 '$1>=0.8 && ($1==0 || ($1-$2)/$1>=0.8){print c,(NR-1)*b-2501,NR*b+2501}' > ${CHR}.${bin}.Zang1817_absbin_ext2k5.bed
  #sed -i 's/-.*\t/0\t/' ${CHR}.${bin}.Zang1817_absbin_ext2k5.bed
  bedtools merge -i ${CHR}.${bin}.Zang1817_absbin_ext2k5.bed |gawk -vOFS="\t" '{$4=$3-$2;if($4>10002){$2=$2+2501;$3=$3-2501;$4=$4-5002;print}}' > ${CHR}.${bin}.Zang1817_absbin_atleasttwobin.bed
done
fi

cat chr??.5k.Zang1817_absbin_atleasttwobin.bed|cut -f4 > Zang1817_absbin_atleasttwobin.length.txt
