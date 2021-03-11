#!/bin/bash
#chr1A chr2A chr3A chr4A chr5A chr6A
dev="/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL"
for CHR in  chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do #chr1A chr2A chr3A chr4A chr5A chr6A
n=$(wc -l < ../GSR_MCL/${CHR}.1M.bed)
#parallel --xapply -j procfile ./bisample_compare_bin_edge_list.sh ../../metadata_cultivar_199.txt ::: $(eval echo {1..$n}) ::: ../../metadata_cultivar_199.txt ::: ${CHR} ::: $(eval echo {1..$n}) ::: homo_undefined_snp_level ::: ~/mapping/fieldergenomecompare/20200424_201_sample_compare
for num in $(eval echo {1..$n});do
while read line;do
echo $line | awk '{if(NF=2) print $0}'| tr ' ' '\n'|cut -d '_' -f2
#echo $line |tr ' ' '\n' | wc -l
done  < ${dev}/mcl_dir/${CHR}.${num}_mcl # | awk '{if ($1>2) print $1}'
cat ../GSR_MCL/mcl_dir/${CHR}.${num}_mcl |  tr '\t' '\n'  | cut -d '_' -f2 | awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($3 in a !=1) print $3}' -  ../../metadata_cultivar_196.txt
#echo $CHR	$num	$a 
#parallel --xapply -j procfile ./bisample_compare_bin_share_num.sh sample_path.txt ::: $(eval echo {1..$n}) ::: ${CHR} ::: $(eval echo {1..$n}) ::: homo_undefined_snp_level ::: ~/mapping/fieldergenomecompare/20200424_201_sample_compare ::: all_sample ::: ../../metadata_cultivar_CNC.txt

done
done
#homo_hmm_snp_level
