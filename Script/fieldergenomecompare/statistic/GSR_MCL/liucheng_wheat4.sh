#!/bin/bash
#for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do #chr1A chr2A chr3A chr4A chr5A chr6A
for CHR in  chr1A chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
n=$(wc -l < ${CHR}.1M.bed)
#parallel --xapply -j procfile ./bisample_compare_bin_edge_list.sh ../../metadata_cultivar_199.txt ::: $(eval echo {1..$n}) ::: ../../metadata_cultivar_199.txt ::: ${CHR} ::: $(eval echo {1..$n}) ::: homo_undefined_snp_level ::: ~/mapping/fieldergenomecompare/20200424_201_sample_compare

parallel --xapply -j procfile_wheat4 ./bisample_compare_bin_edge_list.sh sample_path.txt ::: $(eval echo {1..$n}) ::: ${CHR} ::: $(eval echo {1..$n}) ::: homo_undefined_snp_level.sorted.v1 ::: /data/user/shinyug/HMM_for_yzz_comp/200921/201_sample_compare_HMMdata ::: ../../metadata_cultivar_198.txt
done
#homo_hmm_snp_level
