#!/bin/bash
set -euxo pipefail

# grep 'G[TG]' /data2/rawdata2/sample_metadata/tetra_intro/tetra-introgress.txt|cut -f 2 > WGS_samplelist_l2.txt
# grep 'G[HD]' /data2/rawdata2/sample_metadata/tetra_intro/tetra-introgress.txt|cut -f 2 > D_WGS_samplelist_l2.txt

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
(while read line;do
  gawk 'NR>1{if($1<0.5 || $1>1.3){print 0}else{print 1}}' ~/mapping/fieldergenomecompare/20200416_201_CNV/${line}_2020-04-17162314/${CHR}.1M.norm | /data2/rawdata2/bin/datamash-1.4/datamash transpose
done < metadata_cultivar_samplelist.txt) | /data2/rawdata2/bin/datamash-1.4/datamash transpose > ${CHR}.CNVfilter.txt
done

