#!/usr/bin/env usr/sh

# sum bin to 100k
SM=$1
for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  gawk '{sum+=$4} (NR%100)==0{print sum/100; sum=0}' ../readDepth_ori_1k/$SM/${chr}.1* > $SM/${chr}.1.100k.tmp
  gawk '{sum+=$4} (NR%100)==0{print sum/100; sum=0}' ../readDepth_ori_1k/$SM/${chr}.2* > $SM/${chr}.2.100k.tmp
done

# normlize DP coverage by sample
sum=`gawk '{sum+=$1} END{print sum}' ${SM}/*100k.tmp`
line=`wc -l ${SM}/*100k.tmp|tail -n 1|gawk '{print $1}'`
ave=$(perl -e "print $sum/$line")

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${SM}/${CHR}.1.100k.tmp > ${SM}/${CHR}.100k.norm
  gawk -vOFS="\t" -vave=$ave '{print $1/ave}' ${SM}/${CHR}.2.100k.tmp >> ${SM}/${CHR}.100k.norm
  sed -i '1 s/^.*$/'$SM'/' ${SM}/${CHR}.100k.norm
done

# count ploid bias
while read line;do
  for chr in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
    gawk 'NR==1{print} NR>1{print $1/9*4}' ${line}/${chr}.100k.norm |sponge ${line}/${chr}.100k.norm
  done
done < sample_tetra.txt

# di
for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A;do
  gawk 'NR==1{print} NR>1{print $1/3}' Tu2p0Simu/${chr}.100k.norm |sponge Tu2p0Simu/${chr}.100k.norm
done

# di 
for chr in chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
  gawk 'NR==1{print} NR>1{print $1/3}' AetV4p0Simu/${chr}.100k.norm |sponge AetV4p0Simu/${chr}.100k.norm
done

# join DP by chr
while read line;do paste ../*/${line}.100k.norm > ${line}.100k.join; done < chrlst
mv *.join picture_100k

# join DP info by sample for fingerprint
while read line;do for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do sed '1 s/^.*$/'$chr'/' ${line}/${chr}.1M.norm > ${line}/${chr}.1M.norm.chr done done < sample

while read line;do cd $line; paste *.chr > ${line}.chr;cd ..; done < sample
