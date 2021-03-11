while read line ; do
Chr=`echo $line | awk '{print $1}'`
Pos=`echo $line | awk '{print $2}'`
up=`echo $Pos | awk  '{print $1-300}'`
down=`echo $Pos | awk '{print $1+300}'`

bcftools query -r ${Chr}:${up}-${down} -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' -e 'MAC/N_SAMPLES*2<0.01' /data/user/yangzz/mapping/08.mergeGVCF/196_9_8_dev/${Chr}.bcf.gz
#echo -e "$Chr\t$Pos\t$lis"
done < marker_pos_uniq.txt
