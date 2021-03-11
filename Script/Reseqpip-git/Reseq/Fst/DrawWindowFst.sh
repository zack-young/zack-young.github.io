#!/usr/bin/env sh
#  Guo, Weilong; guoweilong2126.com; 2018-08-03

# --------

PAIR=$1
#PAIR="HA1_vs_LA1"
S1=`echo $PAIR | cut -d"_" -f1`
S2=`echo $PAIR | cut -d"_" -f3`

#size=200 #200k
size=$2

cut -f1 ../../GRP/${S1}.grp > ${S1}.grp
cut -f1 ../../GRP/${S2}.grp > ${S2}.grp

for BCF in ../../10.concatVCF/chr[1-7][ABD].bcf.gz; do
  CHR=`basename $BCF | sed s/.bcf.gz//g`
  echo $CHR
  vcftools --bcf $BCF --weir-fst-pop ${S1}.grp --weir-fst-pop ${S2}.grp  --fst-window-size 200000 --fst-window-step 200000 --out ${CHR}.${PAIR}.${size}k &
  wait_all
done

wait


THR=`gawk 'FNR>1{print $5}' chr[1-7][ABD].DE1_vs_DO1.200k.windowed.weir.fst | st --percentile=95`

# ===
# Draw figure by chromosome

for F in chr[1-7][ABD]*.${PAIR}.${size}k.windowed.weir.fst; do
  PRE=`echo $F | sed s/.windowed.weir.fst//g`
  echo $PRE
  ./DrawWindowFst.R -c $PRE -t $THR &
  sleep 0.2
done

#==
# Draw figure 
for ID in A B D; do
      (cat chr?${ID}.${PAIR}.${size}k.windowed.weir.fst | gawk 'NR==1 || /^chr/' > chr${ID}.${PAIR}.${size}k.windowed.weir.fst
      ./DrawWindowFst_WholeGenome.R -c chr${ID}.${PAIR}.${size}k -w ${size}000 -t $THR ) &
      sleep 1
done
#==

