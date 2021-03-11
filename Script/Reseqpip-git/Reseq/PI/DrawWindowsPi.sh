#!/usr/bin/env sh
#  Guo, Weilong; guoweilong@126.com; 2018-08-03

# --------

PAIR=$1
#PAIR="HA1_vs_LA1"
S1=`echo $PAIR | cut -d"_" -f1`
S2=`echo $PAIR | cut -d"_" -f3`

#size1=20  #20k
#size2=200 #200k

size1=$2
size2=$3


cut -f1 ../../GRP/${S1}.grp > ${S1}.grp
cut -f1 ../../GRP/${S2}.grp > ${S2}.grp

wait

for F in ../../10.concatVCF/chr[1-7][ABD].snp.bcf.gz; do
  echo $F
  CHR=`basename $F | sed s/.bcf.gz//g`
  for Code in ${S1} ${S2}; do
    vcftools --bcf $F --keep ${Code}.grp --window-pi ${size1}000 --out ${CHR}.${Code}.${size1}k &
    vcftools --bcf $F --keep ${Code}.grp --window-pi ${size2}000 --out ${CHR}.${Code}.${size2}k &
  done
  wait_all
done

wait

for id in 1 2 3 4 5 6 7 ;  do
  for code in A B D; do
    CHR=chr${id}${code}
    echo $CHR
    ./DrawWindowsPi.R -c ${CHR}".snp.${S1}" -s ${size1} -S ${size2} &
    ./DrawWindowsPi.R -c ${CHR}".snp.${S2}" -s ${size1} -S ${size2} &
  done 
  wait_all
done

wait


# ----

for id in 1 2 3 4 5 6 7 ;  do
  for code in A B D; do
    CHR=chr${id}${code}
    echo $CHR
    ( gawk -vOFS="\t" 'ARGIND==1{A[$1":"$2]=$5;} ARGIND==2&&(FNR>1){if($1":"$2 in A){print $1, $2, $3, A[$1":"$2], $5, A[$1":"$2]/$5, $5/A[$1":"$2];}}' ${CHR}.snp.${S1}.${size2}k.windowed.pi ${CHR}.snp.${S2}.${size2}k.windowed.pi > ${CHR}.snp.${PAIR}.${size2}k.windowed.pi2; 
    ./DrawWindowsPi_Ratio.R -c ${CHR}".snp.${PAIR}" -s ${size2} ) &
  done
  sleep 1
done

wait

for code in A B D; do
  cat chr[1-7]${code}.snp.${PAIR}.${size2}k.*pi2 > chr${code}.snp.${PAIR}.${size2}k.windowed.pi2
  ./DrawWindowsPi_RatioWholeGenome.R  -c "chr${code}.snp.${PAIR}" -s ${size2}
done

