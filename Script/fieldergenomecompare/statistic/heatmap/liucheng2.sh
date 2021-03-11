#!/bin/bash
num=1
while read lines;do
line_num=`echo $lines|tr ' ' '\n' | wc -l`
#if test ${line_num} != 2 ;then
echo $lines|tr ' ' '\n' | awk -F '\t' 'NR==FNR{a[$1]=$1;next}NR>FNR{if($4"_"$3 in a ==1) print $0}' -  ~/mapping/fieldergenomecompare/metadata_cultivar_196_nongke_10+.txt | cut -f2 > 34Mb_${num}
num=$(( $num +1 ))
#fi
done < ~/mapping/fieldergenomecompare/statistic/GSR_MCL/ALL_mcl_dir/chr2D.35_mcl
#while read line;do
#  sh script_gene.sh -a metadata_all.append.txt -g $line use_samplelist 
#done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/gene.list

#awk -F '\t' 'NR==FNR{a[$3]=$0;next}NR>FNR{if($1 in a ==1) print a[$1]}' metadata_cultivar_198.txt rht4D_sample_list > metadata_cultivar_rht4D.txt

#sh script_region.sh -a metadata_cultivar_196.txt  -c chr2D -s 33452048  -e 34456269 -p ~/mapping/08.mergeGVCF/196_9_8_dev use_samplelist
