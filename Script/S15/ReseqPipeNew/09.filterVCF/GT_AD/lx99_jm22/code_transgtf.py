cat gffcmp.combined.gtf | grep exon | cut -f1,4,5,9 | cut -f1 -d";" |  awk '{print $1, $2, $3, $5}' | sed -e 's/ /\t/g' | sed -e 's/\"//g' > gffcmp.combined.exon.bed 


convert2bed --input=gtf --output=bed < Triticum_aestivum.IWGSC.gtf > ./IWGSC_V1P1_gtf.bed 
