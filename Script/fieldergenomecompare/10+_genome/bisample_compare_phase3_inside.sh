#!/bin/bash
#set -euxo pipefail
line=$1
meta=$2
SAMPLE1=`echo $line|awk -F ',' '{print $1}'`
SAMPLE2=`echo $line|awk -F ',' '{print $2}'`
declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
done < ${meta}

for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do   
    #./muti_func_snp_compare.py --two_diff_level on -i ${SAMPLE1}_${SAMPLE2}/${chr}.1M.delcnv_density --sample1 ${SAMPLE1} --sample2 ${SAMPLE2} 
    #cat ${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_undefined_level ${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV > ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level
    ./ptf_maker.py -p "/data3/user3/wangwx/projs/HMM_for_yzz_comp/201212/201_sample_compare_HMMdata/${SAMPLE1}_${SAMPLE2}" -s ".homo_undefined_snp_level.sorted.v4"> plotfile/plotfile_${SAMPLE1}_${SAMPLE2}
    Rscript ~/R/graph_diff_sim.R -d "~/mapping/fieldergenomecompare/10+_genome" --sample1 ${SAMPLE1} --sample1_name ${arr_3[${SAMPLE1}]} --sample2 ${SAMPLE2} --sample2_name ${arr_3[${SAMPLE2}]} -s ""
done
         #wait
         #./ptf_maker.py -p "${DEV_PATH}/${SAMPLE1}_${SAMPLE2}" -s ".homo_undefined_snp_level"> plotfile/plotfile_${SAMPLE1}_${SAMPLE2}
         #Rscript ~/R/graph_diff_sim.R -d "~/mapping/fieldergenomecompare/202009_11_yaoyy" --sample1 ${SAMPLE1} --sample1_name ${arr_3[${SAMPLE1}]} --sample2 ${SAMPLE2} --sample2_name ${arr_3[${SAMPLE2}]} -s "" 
         #sleep 0.00000000001s
#wait_all
