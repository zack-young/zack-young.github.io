#!/usr/bin/env sh
set -euxo pipefail

S1=BR1
size=1000
cut -f1 /data2/rawdata2/GRP/brittlePheno/${S1}.grp > ${S1}.grp
genefile=/data2/rawdata2/Fst/scripts/DEDO_genelist.txt

if true;then
for BCF in /data/user/yangzz/mapping/08.mergeGVCF/201_final/chr[1-7][ABD].bcf.gz; do
  CHR=`basename $BCF | sed s/.bcf.gz//g`
  #vcftools --bcf $BCF --TajimaD ${size}000 --out ${S1}.${CHR}.${size}k `cat ${S1}.grp | gawk '{printf(" --indv %s", $1)}'` &
  vcftools --bcf $BCF --TajimaD ${size}000 --out ${S1}.${CHR}.${size}k & #`cat ${S1}.grp | gawk '{printf(" --indv %s", $1)}'` &
  wait_all
done
wait

#gawk -vw=${size}000 -vOFS="\t" 'ARGIND==1{a[$1]=$2} ARGIND==2{if(FNR==1){print "CHROM","START","END","TajimaD"}else{if($2+w<a[$1]){e=$2+w}else{e=a[$1]};print $1,$2,e,$4}}' ../chrlength.txt <(cat BR1.chr?A.200k.Tajima.D) | gawk 'NR==1 || /^chr/'> chrA.BR1.200k.Tajima.D.txt
#gawk -vw=${size}000 -vOFS="\t" 'ARGIND==1{a[$1]=$2} ARGIND==2{if(FNR==1){print "CHROM","START","END","TajimaD"}else{if($2+w<a[$1]){e=$2+w}else{e=a[$1]};print $1,$2,e,$4}}' ../chrlength.txt <(cat BR1.chr?B.200k.Tajima.D) | gawk 'NR==1 || /^chr/'> chrB.BR1.200k.Tajima.D.txt
#gawk -vw=${size}000 -vOFS="\t" 'ARGIND==1{a[$1]=$2} ARGIND==2{if(FNR==1){print "CHROM","START","END","TajimaD"}else{if($2+w<a[$1]){e=$2+w}else{e=a[$1]};print $1,$2,e,$4}}' ../chrlength.txt <(cat BR1.chr?D.200k.Tajima.D) | gawk 'NR==1 || /^chr/'> chrD.BR1.200k.Tajima.D.txt
fi
if false;then
THRu=`gawk 'FNR>1 && $4>0{print $4}' chr[AB].BR1.200k.Tajima.D.txt | st --percentile=97.5`
THRl=`gawk 'FNR>1 && $4<0{print $4}' chr[AB].BR1.200k.Tajima.D.txt | st --percentile=2.5`
ID=A
Rscript /data2/rawdata2/Fst/scripts/DrawWindowFst_WholeGenome.R -i chr${ID}.BR1.200k.Tajima.D.txt --threshold="${THRu},${THRl}" -c 1 -g ${genefile} -n 4 -o chr${ID}.BR1.200k.Tajima.D.pdf -T "chr${ID}.BR1.200k.Tajima.D" -a TajimaD --ylolimP="-5" -Y 6
ID=B
Rscript /data2/rawdata2/Fst/scripts/DrawWindowFst_WholeGenome.R -i chr${ID}.BR1.200k.Tajima.D.txt --threshold="${THRu},${THRl}" -c 2 -g ${genefile} -n 4 -o chr${ID}.BR1.200k.Tajima.D.pdf -T "chr${ID}.BR1.200k.Tajima.D" -a TajimaD --ylolimP="-5" -Y 6
THRu=`gawk 'FNR>1 && $4>0{print $4}' chrD.BR1.200k.Tajima.D.txt | st --percentile=97.5`
THRl=`gawk 'FNR>1 && $4>0{print $4}' chrD.BR1.200k.Tajima.D.txt | st --percentile=2.5`
ID=D
Rscript /data2/rawdata2/Fst/scripts/DrawWindowFst_WholeGenome.R -i chr${ID}.BR1.200k.Tajima.D.txt --threshold="${THRu},${THRl}" -c 3 -g ${genefile} -n 4 -o chr${ID}.BR1.200k.Tajima.D.pdf -T "chr${ID}.BR1.200k.Tajima.D" -a TajimaD --ylolimP="-5" -Y 6

qpdf --empty --pages chrA.BR1.200k.Tajima.D.pdf chrB.BR1.200k.Tajima.D.pdf chrD.BR1.200k.Tajima.D.pdf -- BR1.200k.Tajima.D.pdf
fi
