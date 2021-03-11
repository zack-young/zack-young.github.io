while read i;do
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t[%GT\t%AD\t%GQ\t]\n' chr7D.ann.bcf.gz|grep ${i} - >> chen.result
done < chen.txt
