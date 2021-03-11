#!/bin/bash
#./mofaso_third.py -b 1000000 -i /data/user/yangzz/mapping/s1_4/S1/chr1A.raw.vcf.filter  -s 5 -o tes11
#for VCF in /data/user/yangzz/mapping/5xresult/chr??.raw.bcf.filter; do
#  CHR=`basename $VCF | sed s/.raw.bcf.filter//g`
#  echo "${VCF} ${CHR}"
#  ./fouregion.py -i ${CHR} -o ${CHR}.four_filter
#done
for VCF in /data/user/yangzz/mapping/5xresult/calfour/chr??.four_filter; do
  CHR=`basename $VCF | sed s/.four_filter//g`
  echo "${VCF} ${CHR}"
  ./count.py -b 1000000 -i ${VCF} -1 5 -2 7 -s 6  >> total_CNV.count
done

#for VCF in /data/user/yangzz/mapping/5xresult/chr??.raw.bcf.10filter; do
#  CHR=`echo $VCF | sed s/.raw.bcf.10filter//g`
#  echo "${VCF} ${CHR}"
#  ./fivekinds.py -b 1000000 -i ${VCF} -1 5 -2 7 -s 6 -o ${CHR}_10score_CNV &
#done
# cat genotype_count.txt |awk 'ORS=NR%6?"\t":"\n"{print}' |sed 's/|/ /' > genotype_count_line.txt # turn row to line
#echo "total snp" > chr2A.200_600M.10DP_GQ.txt
#less chr2A.200_600M.10DP_GQ | wc -l >> chr2A.200_600M.10DP_GQ.txt
#echo "without ./." >> chr2A.200_600M.10DP_GQ.txt
#less chr2A.200_600M.10DP_GQ|gawk -vOFS="\t" '$5 != "./."&&$6 != "./."&&$7 != "./."'|wc -l >> chr2A.200_600M.10DP_GQ.txt
#echo "5181 lx987" >> chr2A.200_600M.10DP_GQ.txt
#less chr2A.200_600M.10DP_GQ|gawk -vOFS="\t" '$5 != "./."&&$6 != "./."&&$7 != "./."' |gawk -vOFS="\t" '$5 == $6 && $6 != $7 && $6 !={print}' |wc -l >>chr2A.200_600M.10DP_GQ.txt
#echo "5181 3097" >> chr2A.200_600M.10DP_GQ.txt
#less chr2A.200_600M.10DP_GQ|gawk -vOFS="\t" '$5 != "./."&&$6 != "./."&&$7 != "./."' |gawk -vOFS="\t" '$5 != $6 && $6 == $7{print}' |wc -l >>chr2A.200_600M.10DP_GQ.txt
#echo "lx987 5181 3097" >> chr2A.200_600M.10DP_GQ.txt
#less chr2A.200_600M.10DP_GQ|gawk -vOFS="\t" '$5 != "./."&&$6 != "./."&&$7 != "./."' |gawk -vOFS="\t" '$5 == $6 && $6 == $7{print}' |wc -l >>chr2A.200_600M.10DP_GQ.txt
#echo "without similarity" >> chr2A.200_600M.10DP_GQ.txt
#less chr2A.200_600M.10DP_GQ|gawk -vOFS="\t" '$5 != "./."&&$6 != "./."&&$7 != "./."' |gawk -vOFS="\t" '$5 != $6 && $6 != $7 && $5 != $7{print}' |wc -l >>chr2A.200_600M.10DP_GQ.txt

#less chr2A.200_600M.10DP_GQ|gawk -vOFS="\t" '$5 != "./."&&$6 != "./."&&$7 != "./."' |gawk -vOFS="\t" '$5 != $6 && $6 == $7{print}'
