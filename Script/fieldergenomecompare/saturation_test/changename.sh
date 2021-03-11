#!/usr/bin/env sh
# Guo, Weilong; guoweilong@126.com; 2017-10-24
#echo $CHR
chrlst=chr4A.1,chr4A.2,chr4B.1,chr4B.2,chr4D.1,chr4D.2,chr5A.1,chr5A.2,chr5B.1,chr5B.2,chr5D.1,chr5D.2,chr6A.1,chr6A.2,chr6B.1,chr6B.2,chr6D.1,chr6D.2,chr7A.1,chr7A.2,chr7B.1,chr7B.2,chr7D.1,chr7D.2
CHR=`echo ${chrlst} | tr "," " "`
for chr in ${CHR};do
(for DIR in 05_sample 10_sample 25_sample 50_sample 75_sample 85_sample;do
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/picard.jar RenameSampleInVcf \
  INPUT=${DIR}/${chr}.g.vcf.gz \
  OUTPUT=${DIR}/${chr}.g.ch.vcf.gz \
  NEW_SAMPLE_NAME=${DIR}
bcftools index -f -t ${DIR}/${chr}.g.ch.vcf.gz -o ${DIR}/${chr}.g.ch.vcf.gz.tbi
done) &
wait_all
done
