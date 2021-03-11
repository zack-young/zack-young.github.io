#!/usr/bin/env bash
set -euxo pipefail

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
  plink --bcf ../190823/${CHR}.trimed.bcf.gz --make-bed --allow-extra-chr --out ${CHR}_setname --set-missing-var-ids @:# &
  wait_all
done
wait

for i in {1..7};do
for j in A B D;do
  echo chr${i}${j}_setname >> AABBDD.txt
done
done
for i in {1..7};do
for j in A B;do
  echo chr${i}${j}_setname >> AABB.txt
done
done
for i in {1..7};do
for j in D;do
  echo chr${i}${j}_setname >> DD.txt
done
done
for i in AABB DD AABBDD;do
  plink --merge-list ${i}.txt --make-bed --out ${i}_snp --allow-extra-chr
done
#
wait


