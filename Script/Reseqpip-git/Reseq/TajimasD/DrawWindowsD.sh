#!/usr/bin/env sh
# Guo, Weilong; guoweilong2126.com; 2018-08-09


PAIR=$1
#PAIR="HA1_and_LA1"
S1=`echo $PAIR | cut -d"_" -f1`
S2=`echo $PAIR | cut -d"_" -f3`

#size1=20  #20k
#size2=200 #200k

size1=$2
size2=$3


cut -f1 ../../GRP/${S1}.grp > ${S1}.grp
cut -f1 ../../GRP/${S2}.grp > ${S2}.grp


for BCF in ../../10.concatVCF/chr[1-7][ABD].snp.bcf.gz; do
  echo $BCF
  CHR=`basename $BCF | sed s/.snp.bcf.gz//g`
  echo $CHR
  for size in $size1 $size2; do
    echo $CHR; echo $size;
    vcftools --bcf $BCF --TajimaD ${size}000 --out ${S1}.${CHR}.${size}k `cat ${S1}.grp | gawk '{printf(" --indv %s", $1)}'` &
    vcftools --bcf $BCF --TajimaD ${size}000 --out ${S2}.${CHR}.${size}k `cat ${S2}.grp | gawk '{printf(" --indv %s", $1)}'` &
    vcftools --bcf $BCF --TajimaD ${size}000 --out ${S1}_and_${S2}.${CHR}.${size}k `cat ${S1}.grp ${S2}.grp | gawk '{printf(" --indv %s", $1)}'` &
    wait_all
  done
done

wait

for F in `ls *.${size1}k.Tajima.D `;  do
    CHR=`echo $F | cut -d'.' -f1,2`
    echo $CHR;
    ./DrawWindowsD.R -c $CHR -s $size1 -S $size2 & 
    sleep 1
done

#==
# Draw figure
for ID in A B D; do
  for S in ${PAIR} ${S1} ${S2}; do
      (cat ${S}.chr?${ID}.${size2}k.Tajima.D | gawk 'NR==1 || /^chr/' > ${S}.chr${ID}.${size2}k.Tajima.D
      ./DrawWindowsD_WholeGenome.R -c ${S}.chr${ID}.${size2}k -w ${size2}000 ) &
      sleep 1
  done
done

