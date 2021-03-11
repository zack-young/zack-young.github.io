#!/bin/bash
set -euxo pipefail
if false;then
for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do #chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#CHR="chr1A"
n=$(wc -l < ~/mapping/fieldergenomecompare/statistic/GSR_MCL/${CHR}.1M.bed)
parallel --xapply -j procfile --results ALL_sample ./bisample_compare_bin_edge_list.sh sample_path.txt ::: $(eval echo {1..$n}) ::: ${CHR} ::: $(eval echo {1..$n}) ::: homo_undefined_snp_level.sorted.v4 ::: /data3/user3/wangwx/projs/HMM_for_yzz_comp/201212/201_sample_compare_HMMdata ::: ALL_sample ::: ../../metadata_cultivar_196_nongke_10+.txt

done
wait
fi


for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#for CHR in  chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D ;do
n=$(wc -l < ~/mapping/fieldergenomecompare/statistic/GSR_MCL/${CHR}.1M.bed)
parallel --xapply -j procfile --results parallel_result ./bisample_compare_hap_count.sh ::: ${CHR} ::: $(eval echo {1..$n}) ::: CNL
#if true;then
done
wait 
#for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#sort -nk2,2 ${CHR}_all_hap  | sponge ${CHR}_all_hap
#done
