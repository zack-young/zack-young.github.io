#!/bin/bash
set -x

while getopts 'c:v:s:m:d:f:n:' opt
do
    case "$opt" in
    c) SAMPLE1=$OPTARG ;;
    s) SAMPLE2=$OPTARG ;;
    esac
done

if [[ -z "${SAMPLE1}" ]]; then
    echo "No first sample provided"; exit 1
fi
if [[ -z "${SAMPLE2}" ]]; then
    echo "No second sample provided."; exit 1
fi
if [[ ! -d "tmp" ]]; then
     mkdir tmp
fi
if [[ ! -d "tmp/${SAMPLE1}" ]]; then
     mkdir tmp/${SAMPLE1}
fi
if [[ ! -d "tmp/${SAMPLE2}" ]]; then
     mkdir tmp/${SAMPLE2}
fi
if [[ ! -d "tmp/${SAMPLE1}_${SAMPLE2}" ]]; then
     mkdir tmp/${SAMPLE1}_${SAMPLE2}
fi



for i in {${SAMPLE1},${SAMPLE2}}; do
              mkdir CNV_heatmap/${i}
              echo ${i} > /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt
              while read i;do
                while read j;do
                  sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/script_ref2.sh $i $j &
                  wait_all
                done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/chrlst
              done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt
              # norm       
              while read i;do
                sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/normDPbySample.sh $i "1k.DP" "norm" 4 &
                wait_all
              done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt
              wait
              while read line;do
                sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/script_1M.sh $line
              done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt
done
wait
for num in {1..7}; do
    for i in {A,B,D}; do
        #1
        for sample in {${SAMPLE1},${SAMPLE2}}; do
            if [[ ! -s "tmp/${SAMPLE1}/chr${num}${i}_snp_${item}" ]];then
                 if [[ -d "CNV_heatmap/${sample}" ]]&&[[ -s "CNV_heatmap/${sample}/chr${num}${i}.1M.norm" ]]; then
                      ./muti_func_snp_compare.py -b 1000000 -i CNV_heatmap/${sample}/chr${num}${i}.1M.norm --mask_cnv on -o tmp/${sample}/chr${num}${i}.mask_CNV &
                 else
                      ./muti_func_snp_compare.py -b 1000000 -i /data2/rawdata2/readDepth/${sample}/chr${num}${i}.1M.norm --mask_cnv on -o tmp/${sample}/chr${num}${i}.mask_CNV &
                 fi
            fi
        done
        wait
        #2
        awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE2"_CNV"}' SAMPLE2="$SAMPLE2" tmp/${SAMPLE1}/chr${num}${i}.mask_CNV tmp/${SAMPLE2}/chr${num}${i}.m    ask_CNV > tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE2}to${SAMPLE1}_own_CNV
        awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE1"_CNV"}' SAMPLE1="$SAMPLE1" tmp/${SAMPLE2}/chr${num}${i}.mask_CNV tmp/${SAMPLE1}/chr${num}${i}.m    ask_CNV > tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV
        awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==1) print $1"\t"$2"\t"$3"\t""both_CNV"}' tmp/${SAMPLE1}/chr${num}${i}.mask_CNV tmp/${SAMPLE2}/chr${num}${i}.mask_CNV > tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}both_CNV
        #3
        cat tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE2}to${SAMPLE1}_own_CNV tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}both_CNV | sort -u -nk2,2 > tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV
        wait_all
    done
done

