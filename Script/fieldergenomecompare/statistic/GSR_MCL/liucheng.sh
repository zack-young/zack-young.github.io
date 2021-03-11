#!/bin/bash
for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do #chr1A chr2A chr3A chr4A chr5A chr6A
n=$(wc -l < ${CHR}.1M.bed)
parallel --xapply -j procfile ./bisample_compare_hap_length_type_density.sh ::: ${CHR} ::: $(eval echo {1..$n}) 

#parallel --xapply -j procfile ./bisample_compare_bin_edge_list.sh sample_path_C.txt ::: $(eval echo {1..$n}) ::: ${CHR} ::: $(eval echo {1..$n}) ::: homo_undefined_snp_level ::: ~/mapping/fieldergenomecompare/20200424_201_sample_compare ::: CNC_sample ::: ../../metadata_cultivar_CNC.txt
done
#homo_hmm_snp_level
#parallel --xapply -j procfile ./bisample_compare_hap_length_type_density.sh ::: chr1A ::: 1 ::: $(eval cat sample_list_CNC.txt)
