#!/bin/bash
#chr1A chr2A chr3A chr4A chr5A chr6A
dev="/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL"
parallel -j procfile sh bisample_compare_consec_region.sh ::: $(eval cat sample_path_noNK.txt) ::: /data3/user3/wangwx/projs/HMM_for_yzz_comp/GaussianHMM_pipeline/201_sample_compare_HMMdata ::: ../../metadata_cultivar_196_nongke_10+.txt 

#homo_hmm_snp_level
