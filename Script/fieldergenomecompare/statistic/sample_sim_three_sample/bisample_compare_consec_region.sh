#!/bin/bash
#set -euxo pipefail
sp_line=$1
dev=$2
meta_data=$3
out_dir=$4


declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
done < ${meta_data}
#SAMPLE_list=${SAMPLE1}
#arr=(`echo ${SAMPLE_list} | tr ',' ' '`)
#arr_3=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
#declare -A arr_2
#while read line;do
#arr_2[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $1}'`
#done < ${meta_data}

num_tmp=0
line_num=1
num_diff=0
dev="/data3/user3/wangwx/projs/HMM_for_yzz_comp/GaussianHMM_pipeline/201_sample_compare_HMMdata"
#dir="/data/user/yangzz/mapping/fieldergenomecompare/statistic/Rht8"
if true;then
           #myjob=`jobs|wc -l`
           SAMPLE1=`echo $sp_line | awk -F ',' '{print $1}'`
           SAMPLE2=`echo $sp_line | awk -F ',' '{print $2}'`
#awk -F '_' '{if(NF >= 2) { for(i=2;i<=NF-1;i++) printf $i"\t";printf $NF"\n"}}'
           #num_low=0
           #num_diff=0
           #num_cnv=10
           #echo $num_tmp

           #num1=`cat ${dev}/${SAMPLE1}_${SAMPLE2}/chr*.homo_undefined_snp_level | grep -E 'low|both' | sort -k 1,1 -k 2,2n| bedtools merge -i - -c 4 -o distinct|awk '{n=($3-$2)/1000000} {print n}'|awk -F '.' '{if ($2>0) {n=$1+1} else{n=$1} {print n}}'|awk 'BEGIN{sum=0}{sum+=$1}END{print sum/NR}'`
           row=`cat /data3/user3/wangwx/projs/HMM_for_yzz_comp/GaussianHMM_pipeline/201_sample_compare_HMMdata/${SAMPLE1}_${SAMPLE2}/chr*homo_undefined_snp_level.sorted.v1 | grep -E 'low|both'|wc -l`
           if [[ "$row" > 0 ]];then
           num2=`cat ${dev}/${SAMPLE1}_${SAMPLE2}/chr*.homo_undefined_snp_level | grep -E 'low|both' |wc -l|awk '{print $1/14075}'`
           num1=`cat /data3/user3/wangwx/projs/HMM_for_yzz_comp/GaussianHMM_pipeline/201_sample_compare_HMMdata/${SAMPLE1}_${SAMPLE2}/chr*homo_undefined_snp_level.sorted.v1 | grep -E 'low|both' | sort -k 1,1 -k 2,2n| bedtools merge -i - -c 4 -o distinct|awk '{n=($3-$2)/1000000} {print n}'|awk -F  '.' '{if ($2>0) {n=$1+1} else{n=$1} {print n}}'|grep -vw 1  |awk 'BEGIN{sum=0} {sum+=$1}END{if (NR!=0) {print sum/NR}}'`

           #cat /data3/user3/wangwx/projs/HMM_for_yzz_comp/GaussianHMM_pipeline/201_sample_compare_HMMdata/YM03_YM02/chr*homo_undefined_snp_level.sorted.v1 | grep -E 'low|both' | sort -k 1,1 -k 2,2n| bedtools merge -i - -c 4 -o distinct|awk '{n=($3-$2)/1000000} {print n}'|awk -F  '.' '{if ($2>0) {n=$1+1} else{n=$1} {print n}}'|sort -n  |awk 'BEGIN{sum=0} {sum+=$1 ;if (sum > 4086) {print $1}}'



           if [[ "$num1" > 0 ]]; then 
           echo -e "${arr_3[${SAMPLE1}]}\t${arr_3[${SAMPLE2}]}\t$num1\t$num2" >> simi_binlen.txt
           fi
           fi
           #if [ "$num_diff" = 0 ]&&[ "$num_low" = 0 ] ;then
           #  final_num=$num_cnv
           #else
           #  final_num=`awk 'BEGIN{printf "%.2f\n",'$num_cnv'+'$num_low'/('$num_low'+'$num_diff')}'`
           #fi
           #if [ "$final_num" != '10.00' ]; then
           #  echo -e "${arr_3[${SAMPLE1}]}\t${arr_3[${SAMPLE2}]}\t$final_num" >> ${out_dir}/${chr}_${suffix2}
           #fi
           #if [[ $myjob > 70 ]];then
           #   wait_all
           #fi
line_num=$(( $line_num + 1 ))
fi
wait
#paste `awk '{print "'$dir'""/""'$chr'""_"$1"_""'$suffix2'""_dist"}'  ../metadata_cultivar_final.txt |tr '\n' ' '|sed s/,$//g` > ${dir}/${chr}_combine_${suffix2}_dist
##
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#cat  dev/${chr}_*_homo_level |grep 'low'> ${chr}_combine_homo_low
#done
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_deletion_split"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_deletion_split
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
        #num=${chr:3:2}
        #./count_CNV.py --count on --pop_num 7021 -i ${chr}_combine_homo_low |sort -nk1,1 >  ${chr}_combine_homo_low_count
        #./count_CNV.py --compensent on --chrom ${num} -i ${chr}_combine_homo_low_count -o ${chr}.combine_homo_low_compensent
#done

