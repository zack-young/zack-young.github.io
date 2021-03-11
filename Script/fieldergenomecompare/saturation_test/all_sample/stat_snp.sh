#!/usr/bin/env sh
# Guo, Weilong; guoweilong@126.com; 2017-10-24
#echo $CHR
for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
(for DIR in 05_sample 10_sample 25_sample 50_sample 75_sample 85_sample s30;do
homonum=`bcftools view -v snps ${chr}.bcf.gz| bcftools query -s ${DIR} -f "[%GT\t]\n"| grep '1/1'|wc -l`
hetnum=`bcftools view -v snps ${chr}.bcf.gz | bcftools query -s ${DIR} -f "[%GT\t]\n" | grep '0/1'|wc -l`
echo ${homonum} ${hetnum} ${DIR} >> ${chr}.stat
done) &
done
