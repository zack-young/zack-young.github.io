#!/usr/bin/env sh
# Guo, Weilong; guoweilong@126.com; 2017-10-24
#echo $CHR
for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
for DIR in s3097-1 LX987 S140;do
#a=`grep 'total' ${DIR}/chr*.2.*flagstat | cut -d ' ' -f1|cut -d ':' -f2|awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'`
#b=`grep 'total' ${DIR}/chr*.1.*flagstat | cut -d ' ' -f1|cut -d ':' -f2|awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'`
#echo $a $b|awk '{print ($1+$2)*150/14547261565}' #14547261565
#(grep -v 'SRR512' ${chr}.stat | sponge ${chr}.stat
#DIR=/data3/user3/publicdata/ReseqData/workflowFile/SRR512/06.callsnp
(homonum1=`bcftools view -v snps -s $DIR ${chr}.bcf.gz| bcftools query  -f "[%GT\t]\n"| grep '1/1'|wc -l`
hetnum1=`bcftools view -v snps -s $DIR ${chr}.bcf.gz | bcftools query  -f "[%GT\t]\n" | grep '0/1'|wc -l`
echo  ${homonum1} ${hetnum1} ${DIR} >> ${chr}.stat
) &
done
done
