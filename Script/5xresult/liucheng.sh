#!/bin/bash

for vcf in *raw.bcf; do
  ../secondstep ${vcf}
done
#./secondstep chr1A.1.raw.vcf 
#./secondstep chr7D.2.raw.vcf


#for vcf in *.2.raw.vcf.filter; do
#  awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' ${vcf} > ${vcf}.4b
#done
#awk '{print $1"\t"$2+400000000"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' chr1A.2.raw.vcf.filter > chr1A.2.vcf.GT.filter.4b


#for VCF in /data/user/yangzz/mapping/deepth/chr*.1.raw.vcf.filter; do
#  CHR=`basename $VCF | sed s/.raw.vcf.filter//g | uniq`
#  cat ${VCF} ${CHR}.2.raw.vcf.filter.4b > ${CHR}.vcf.GT.filter
#done


#rm /data/user/yangzz/mapping/deepth/rawdata/???M.graphy
#for VCF in /data/user/yangzz/mapping/5xresult/chr??.find.filter; do
#  CHR=`echo $VCF | sed s/.find.filter//g`
#  awk '{print $5}' ${VCF} > ${CHR}_five
#done

#for VCF in /data/user/yangzz/mapping/5xresult/5181/chr*.raw.bcf.filter; do
#  CHR=`echo $VCF | sed s/.raw.bcf.filter//g | uniq`
#  ./density.py -b 1000000 -i ${VCF} -s 5 -o ${CHR}.find.filter &
#done
#
#for VCF in /data/user/yangzz/mapping/5xresult/lx987/chr*.raw.bcf.filter; do
#  CHR=`echo $VCF | sed s/.raw.bcf.filter//g | uniq`
#  ./density.py -b 1000000 -i ${VCF} -s 5 -o ${CHR}.find.filter &
#done
#
#for VCF in /data/user/yangzz/mapping/5xresult/3097/chr*.raw.bcf.filter; do
#  CHR=`echo $VCF | sed s/.raw.bcf.filter//g | uniq`
#  ./density.py -b 1000000 -i ${VCF} -s 5 -o ${CHR}.find.filter &
#done

#for VCF in /data/user/yangzz/mapping/5xresult/chr*.vcf.GT.filter; do
#  CHR=`echo $VCF | sed s/.vcf.GT.filter//g | uniq`
#  ./fivekinds.py -b 1000000 -i ${VCF}  -1 5 -2 7 -s 6 -o ${CHR}.find.filter &
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
#./lastgeno.py -b 5000000 -i chr7D.vcf.GT.filter -m 3981 -l 251 -d 0.2 -o chr7D.find.filter

#for VCF in /data/user/yangzz/mapping/5xresult/chr*.vcf.GT.filter; do
#  CHR=`echo $VCF | sed s/.vcf.GT.filter//g | uniq`
#  ./threelevel.py -b 1000000 -i ${VCF} -m 708 -l 50  -1 5 -2 7 -s 6 -o ${CHR}.find.filter &
#done


#for VCF in /data/user/yangzz/mapping/5xresult/lx987/chr*.raw.bcf.filter; do
#  CHR=`echo $VCF | sed s/.raw.bcf.filter//g | uniq`
#  ./one_density.py -b 1000000 -i ${VCF} -m 8400 -l 4000 -s 5 -o ${CHR} &
#done
#
#for VCF in /data/user/yangzz/mapping/5xresult/3097/chr*.raw.bcf.filter; do
#  CHR=`echo $VCF | sed s/.raw.bcf.filter//g | uniq`
#  ./one_density.py -b 1000000 -i ${VCF} -m 8700 -l 3900 -s 5 -o ${CHR} &
#done
#
#for VCF in /data/user/yangzz/mapping/5xresult/5181/chr*.raw.bcf.filter; do
#  CHR=`echo $VCF | sed s/.raw.bcf.filter//g | uniq`
#  ./one_density.py -b 1000000 -i ${VCF} -m 8400 -l 3500 -s 5 -o ${CHR} &
#done


#for VCF in /data/user/yangzz/mapping/5xresult/chr*.raw.bcf.filter; do
#  CHR=`echo $VCF | sed s/.raw.bcf.filter//g | uniq`
#  ./mofaso.py -b 1000000 -i ${VCF} -l 5000 -m 12000 -1 5 -2 7 -s 6 -o ${CHR}.find.filter &
#  ./mofaso.py -b 1000000 -i ${VCF} -1 5 -2 7 -s 6 -o ${CHR}.filter &
#done
#
#for VCF in /data/user/yangzz/mapping/5xresult/chr*.raw.bcf.filter; do
#  CHR=`echo $VCF | sed s/.raw.bcf.filter//g | uniq`
#  ./mofaso_third.py -b 1000000 -i ${VCF} -l 5000 -m 12000 -1 5 -2 7 -s 6 -o ${CHR}.find.filter &
#  ./mofaso_define.py -b 1000000 -i ${VCF} -1 5 -2 7 -s 6 -o ${CHR}.filter1 &
#done


#for VCF in /data/user/yangzz/mapping/5xresult/chr??.vcf.GT.filter; do
#  CHR=`echo ${VCF} | sed s/.vcf.GT.filter//g`
#  ./levelandhete.py -b 1000000 -i ${VCF} -1 5 -2 7 -s 6 -o ${CHR}.5x1M &
#  ./levelandhete.py -b 5000000 -i ${VCF} -1 5 -2 7 -s 6 -o ${CHR}.5x5M &
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
#done

#for VCF in chr??.vcf.GT.filter; do
#  CHR=`echo ${VCF} | sed s/.vcf.GT.filter//g`
#  ./levelandhete.py -b 5000000 -i ${VCF} -m 3981 -l 251 -1 5 -2 7 -s 6 -o ${CHR}.5x5M
#done


#for BCF in /data2/rawdata2/mergeFile/10.concatVCF/chr[6,7]?.snp.bcf.gz; do
#  CHR=`basename ${BCF} | sed s/.snp.bcf.gz//g`
#  bcftools view -s LX987 ${BCF} -Oz -o /data/user/yangzz/mapping/5xresult/lx987/${CHR}.raw.vcf &
#  bcftools view -s s3097-1 ${BCF} -Oz -o /data/user/yangzz/mapping/5xresult/3097/${CHR}.raw.vcf &
#  bcftools view -s S140 ${BCF} -Oz -o /data/user/yangzz/mapping/5xresult/5181/${CHR}.raw.vcf &
#done

#for GRAPHY in /data/user/yangzz/mapping/deepth/rawdata/????.graphy; do
#  DIR=`basename ${GRAPHY} | sed s/.graphy//g`
#  echo ${GRAPHY}
#  ./raw_data_wash.py -i ${GRAPHY} -o /data/user/yangzz/mapping/deepth/rawdata/dir${DIR}/
#done
#./mapping123.py -b 5000000 -i chr1B.vcf.GT.filter -d 0.1 -o chr1B.find.filter
#./mapping123.py -b 5000000 -i chr1D.vcf.GT.filter -d 0.1 -o chr1D.find.filter
#./mapping123.py -b 5000000 -i chr2A.vcf.GT.filter -d 0.1 -o chr2A.find.filter
#./mapping123.py -b 5000000 -i chr2B.vcf.GT.filter -d 0.1 -o chr2B.find.filter
#./mapping123.py -b 5000000 -i chr2D.vcf.GT.filter -d 0.1 -o chr2D.find.filter
#./mapping123.py -b 5000000 -i chr7B.vcf.GT.filter -d 0.1 -o chr7B.find.filter
#./mapping123.py -b 5000000 -i chr7D.vcf.GT.filter -d 0.1 -o chr7D.find.filter
#./countcolor.py -i combine1M -o sta1M
#./countcolor.py -i combine2M -o sta2M
#./countcolor.py -i combine3M -o sta3M
#./countcolor.py -i combine4M -o sta4M
#./countcolor.py -i combine5M -o sta5M
