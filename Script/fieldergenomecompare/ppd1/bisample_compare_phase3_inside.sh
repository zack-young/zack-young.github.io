#!/bin/bash
#set -euxo pipefail
line=$1
chr=$2
SAMPLE1=`echo $line|awk -F ',' '{print $1}'`
SAMPLE2=`echo $line|awk -F ',' '{print $2}'`
    #./muti_func_snp_compare.py --two_diff_level on -i ${SAMPLE1}_${SAMPLE2}/${chr}.1M.density --sample1 ${SAMPLE1} --sample2 ${SAMPLE2} --chromosome ${chr}
    line=`cat ${SAMPLE1}_${SAMPLE2}/${chr}.1.${SAMPLE1}to${SAMPLE2}all_CNV|wc -l`
    if test $line = 0;then
    cat ${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_undefined_level  > ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level
    else
    type=`cut -f4 ${SAMPLE1}_${SAMPLE2}/${chr}.1.${SAMPLE1}to${SAMPLE2}all_CNV`
    awk '{print $1"\t"$2"\t"$3"\t"type}' type="$type" ${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_undefined_level > ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level
    fi
         #wait
         #./ptf_maker.py -p "${DEV_PATH}/${SAMPLE1}_${SAMPLE2}" -s ".homo_undefined_snp_level"> plotfile/plotfile_${SAMPLE1}_${SAMPLE2}
         #Rscript ~/R/graph_diff_sim.R -d "~/mapping/fieldergenomecompare/202009_11_yaoyy" --sample1 ${SAMPLE1} --sample1_name ${arr_3[${SAMPLE1}]} --sample2 ${SAMPLE2} --sample2_name ${arr_3[${SAMPLE2}]} -s "" 
         #sleep 0.00000000001s
#wait_all
