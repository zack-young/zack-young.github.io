#!/bin/bash

#for vcf in *raw.vcf; do
#  ./secondstep ${vcf}
#done
#./secondstep chr7D.2.raw.vcf


#for vcf in *.2.raw.vcf.filter; do
#  awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"}' ${vcf} > ${vcf}.4b
#done
#nohup &
#awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' chr1A.2.raw.vcf.filter > chr1A.2.vcf.GT.filter.4b


#for VCF in /data/user/yangzz/mapping/deepth/chr*.1.raw.vcf.filter; do
#  CHR=`basename $VCF | sed s/.1.raw.vcf.filter//g | uniq`
#  cat ${VCF} ${CHR}.2.raw.vcf.filter.4b > ${CHR}.vcf.GT.filter
#done
#cat chr5B.1.raw.vcf.filter chr5B.2.vcf.GT.filter.4b > chr5B.vcf.GT.filter

#for VCF in /data/user/yangzz/mapping/deepth/chr*.vcf.GT.filter; do
#  CHR=`basename $VCF | sed s/.vcf.GT.filter//g | uniq`
#  ./lastgeno.py -b 3000000 -i ${VCF} -m 90 -l 5 -d 0.00092 -1 14 -2 6 -s 10 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.2x3M.graphy
#  ./lastgeno.py -b 4000000 -i ${VCF} -m 97 -l 8 -d 0.00085 -1 14 -2 6 -s 10 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.2x4M.graphy
#  ./lastgeno.py -b 5000000 -i ${VCF} -m 114 -l 11 -d 0.00065 -1 14 -2 6 -s 10 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.2x5M.graphy
#  ./lastgeno.py -b 2000000 -i ${VCF} -m 142 -l 5 -d 0.00027 -1 15 -2 7 -s 11 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.3x2M.graphy
#  ./lastgeno.py -b 3000000 -i ${VCF} -m 186 -l 8 -d 0.00027 -1 15 -2 7 -s 11 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.3x3M.graphy
#  ./lastgeno.py -b 4000000 -i ${VCF} -m 235 -l 12 -d 0.00062 -1 15 -2 7 -s 11 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.3x4M.graphy
#  ./lastgeno.py -b 5000000 -i ${VCF} -m 235 -l 12 -d 0.2 -1 15 -2 7 -s 11 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.3x5M.graphy
#  ./lastgeno.py -b 2000000 -i ${VCF} -m 296 -l 8 -d 0.2 -1 16 -2 8 -s 12 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.4x2M.graphy
#  ./lastgeno.py -b 3000000 -i ${VCF} -m 413 -l 11 -d 0.00012 -1 16 -2 8 -s 12 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.4x3M.graphy
#  ./lastgeno.py -b 4000000 -i ${VCF} -m 551 -l 15 -d 0.2 -1 16 -2 8 -s 12 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.4x4M.graphy
#  ./lastgeno.py -b 5000000 -i ${VCF} -m 551 -l 15 -d 0.2 -1 16 -2 8 -s 12 -o /data/user/yangzz/mapping/deepth/rawdata/${CHR}.4x5M.graphy
#done
#./lastgeno.py -b 2000000 -i chr1A.vcf.GT.filter -m 1585 -l 100 -d 0.2 -1 5 -2 13 -s 9 -o chr1A.find.filter   #-1 3097 -2 987 -s 5181
#./lastgeno.py -b 2000000 -i chr7D.vcf.GT.filter -m 1585 -l 100 -d 0.2 -o chr7D.find.filter


#rm /data/user/yangzz/mapping/deepth/rawdata/???M.graphy
#for VCF in /data/user/yangzz/mapping/deepth/rawdata/chr*M.graphy; do
#  CHR=`basename $VCF | sed s/chr...//g | uniq`
#  cat ${VCF} | awk '{print $1"\t"$2"\t"$3"\t""5181""\t"$7}' >> /data/user/yangzz/mapping/deepth/rawdata/${CHR}
#done


for VCF in chr??.vcf.GT.filter; do
  CHR=`echo ${VCF} | sed s/.vcf.GT.filter//g`
  ./levelandhete.py -b 5000000 -i ${VCF} -1 5 -2 7 -s 6 -o ${CHR}.5x5M
#  ./levelandhete.py -b 2000000 -i ${VCF} -1 5 -2 13 -s 9 -o ${CHR}.1x2M
#  ./levelandhete.py -b 3000000 -i ${VCF} -1 5 -2 13 -s 9 -o ${CHR}.1x3M
#  ./levelandhete.py -b 4000000 -i ${VCF} -1 5 -2 13 -s 9 -o ${CHR}.1x4M
#  ./levelandhete.py -b 5000000 -i ${VCF} -1 5 -2 13 -s 9 -o ${CHR}.1x5M
#  ./levelandhete.py -b 1000000 -i ${VCF} -1 6 -2 14 -s 10 -o ${CHR}.2x1M
#  ./levelandhete.py -b 2000000 -i ${VCF} -1 6 -2 14 -s 10 -o ${CHR}.2x2M
#  ./levelandhete.py -b 3000000 -i ${VCF} -1 6 -2 14 -s 10 -o ${CHR}.2x3M
#  ./levelandhete.py -b 4000000 -i ${VCF} -1 6 -2 14 -s 10 -o ${CHR}.2x4M
#  ./levelandhete.py -b 5000000 -i ${VCF} -1 6 -2 14 -s 10 -o ${CHR}.2x5M
#  ./levelandhete.py -b 1000000 -i ${VCF} -1 7 -2 15 -s 10 -o ${CHR}.3x1M
#  ./levelandhete.py -b 2000000 -i ${VCF} -1 7 -2 15 -s 10 -o ${CHR}.3x2M
#  ./levelandhete.py -b 3000000 -i ${VCF} -1 7 -2 15 -s 10 -o ${CHR}.3x3M
#  ./levelandhete.py -b 4000000 -i ${VCF} -1 7 -2 15 -s 10 -o ${CHR}.3x4M
#  ./levelandhete.py -b 5000000 -i ${VCF} -1 7 -2 15 -s 10 -o ${CHR}.3x5M
#  ./levelandhete.py -b 1000000 -i ${VCF} -1 8 -2 16 -s 11 -o ${CHR}.4x1M
#  ./levelandhete.py -b 2000000 -i ${VCF} -1 8 -2 16 -s 11 -o ${CHR}.4x2M
#  ./levelandhete.py -b 3000000 -i ${VCF} -1 8 -2 16 -s 11 -o ${CHR}.4x3M
#  ./levelandhete.py -b 4000000 -i ${VCF} -1 8 -2 16 -s 11 -o ${CHR}.4x4M
#  ./levelandhete.py -b 5000000 -i ${VCF} -1 8 -2 16 -s 11 -o ${CHR}.4x5M
done

#./mapping123.py -b 5000000 -i chr1B.vcf.GT.filter -d 0.1 -o chr1B.find.filter
#./mapping123.py -b 5000000 -i chr1D.vcf.GT.filter -d 0.1 -o chr1D.find.filter
#./mapping123.py -b 5000000 -i chr2A.vcf.GT.filter -d 0.1 -o chr2A.find.filter
#./mapping123.py -b 5000000 -i chr2B.vcf.GT.filter -d 0.1 -o chr2B.find.filter
#./mapping123.py -b 5000000 -i chr2D.vcf.GT.filter -d 0.1 -o chr2D.find.filter
#./mapping123.py -b 5000000 -i chr3A.vcf.GT.filter -d 0.1 -o chr3A.find.filter
#./mapping123.py -b 5000000 -i chr3B.vcf.GT.filter -d 0.1 -o chr3B.find.filter
#./mapping123.py -b 5000000 -i chr3D.vcf.GT.filter -d 0.1 -o chr3D.find.filter
#./mapping123.py -b 5000000 -i chr4A.vcf.GT.filter -d 0.1 -o chr4A.find.filter
#./mapping123.py -b 5000000 -i chr4B.vcf.GT.filter -d 0.1 -o chr4B.find.filter
#./mapping123.py -b 5000000 -i chr4D.vcf.GT.filter -d 0.1 -o chr4D.find.filter
#./mapping123.py -b 5000000 -i chr5A.vcf.GT.filter -d 0.1 -o chr5A.find.filter
#./mapping123.py -b 5000000 -i chr5B.vcf.GT.filter -d 0.1 -o chr5B.find.filter
#./mapping123.py -b 5000000 -i chr5D.vcf.GT.filter -d 0.1 -o chr5D.find.filter
#./mapping123.py -b 5000000 -i chr6A.vcf.GT.filter -d 0.1 -o chr6A.find.filter
#./mapping123.py -b 5000000 -i chr6B.vcf.GT.filter -d 0.1 -o chr6B.find.filter
#./mapping123.py -b 5000000 -i chr6D.vcf.GT.filter -d 0.1 -o chr6D.find.filter
#./mapping123.py -b 5000000 -i chr7A.vcf.GT.filter -d 0.1 -o chr7A.find.filter
#./mapping123.py -b 5000000 -i chr7B.vcf.GT.filter -d 0.1 -o chr7B.find.filter
#./mapping123.py -b 5000000 -i chr7D.vcf.GT.filter -d 0.1 -o chr7D.find.filter
#./countcolor.py -i combine1M -o sta1M
#./countcolor.py -i combine2M -o sta2M
#./countcolor.py -i combine3M -o sta3M
#./countcolor.py -i combine4M -o sta4M
#./countcolor.py -i combine5M -o sta5M
