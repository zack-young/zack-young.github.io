#!/bin/bash
#set -euxo pipefail
line=$1
SAMPLE1=`echo $line|awk -F ',' '{print $1}'`
SAMPLE2=`echo $line|awk -F ',' '{print $2}'`
for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
    ./muti_func_snp_compare.py --two_diff_level on -i ${SAMPLE1}_${SAMPLE2}/${chr}.1M.delcnv_density --sample1 ${SAMPLE1} --sample2 ${SAMPLE2} --chromosome ${chr}
    cat ${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_undefined_level ${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV > ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level
done
         #wait
         #./ptf_maker.py -p "${DEV_PATH}/${SAMPLE1}_${SAMPLE2}" -s ".homo_undefined_snp_level"> plotfile/plotfile_${SAMPLE1}_${SAMPLE2}
         #Rscript ~/R/graph_diff_sim.R -d "~/mapping/fieldergenomecompare/202009_11_yaoyy" --sample1 ${SAMPLE1} --sample1_name ${arr_3[${SAMPLE1}]} --sample2 ${SAMPLE2} --sample2_name ${arr_3[${SAMPLE2}]} -s "" 
         #sleep 0.00000000001s
#wait_all
