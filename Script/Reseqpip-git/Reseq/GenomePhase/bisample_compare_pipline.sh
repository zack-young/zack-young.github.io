#!/bin/bash
set -x

while getopts 'c:v:s:m:d:f:n:' opt
do
    case "$opt" in
    c) SAMPLE1=$OPTARG ;;
    v) SAMPLE1_name=$OPTARG ;;
    s) SAMPLE2=$OPTARG ;;
    m) SAMPLE2_name=$OPTARG ;;
    d) total=$OPTARG ;;
    f) first=$OPTARG ;;
    n) second=$OPTARG ;;
    esac
done

if [[ -z "${SAMPLE1}" ]]; then
    echo "No first sample provided"; exit 1
fi
if [[ -z "${SAMPLE2}" ]]; then
    echo "No second sample provided."; exit 1
fi
if [[ -z "${SAMPLE1_name}" ]]; then
    echo "No first sample real name provided"; exit 1
fi
if [[ -z "${SAMPLE2_name}" ]]; then
    echo "No second sample real name provided."; exit 1
fi
if [[ -z "${total}" ]]; then
    echo "No total sample bcf name provided."; exit 1
fi
if [[ -z "${first}" ]]; then
    echo "No first sample column provided."; exit 1
fi
if [[ -z "${second}" ]]; then
    echo "No second sample column provided."; exit 1
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
if [[ ! -d "${total}" ]]; then
     mkdir ${total}
fi
if [[ ! -d "tmp/${SAMPLE1}_${SAMPLE2}" ]]; then
     mkdir tmp/${SAMPLE1}_${SAMPLE2}
fi

#for VCF in /data/user/yangzz/mapping/S15/ReseqPipeNew/08.mergeGVCF/${total}/chr??.bcf.gz; do
#    CHR=`basename $VCF | sed s/.bcf.gz//g`
#    echo "${VCF}"
#    #bcftools view  -v indels -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > GT_AD/${CHR}_indel_akcp.gt &
#    if [[ ! -s "${total}/${CHR}_snp.gt" ]]; then
#         bcftools view  -v snps -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > ${total}/${CHR}_snp.gt &
#    fi
#    if [[ ! -s "${total}/${CHR}_indel.gt" ]]; then
#         bcftools view  -v indels -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > ${total}/${CHR}_indel.gt &
#    fi
#    wait_all
#done
#wait
#for i in {${SAMPLE1},${SAMPLE2}}; do
#    if [[ ! -d "/data2/rawdata2/readDepth/${i}" ]]; then 
#         if [[ ! -d "CNV_heatmap/${i}" ]]||[[ ! -n "`find /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap -maxdepth 1 -name 'chr*'`" ]]; then
#              mkdir CNV_heatmap/${i}
#              echo ${i} > /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt
#              while read i;do
#                while read j;do
#                  sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/script.sh $i $j &
#                  wait_all  
#                done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/chrlst 
#              done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt 
#              # norm       
#              while read i;do 
#                sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/normDPbySample.sh $i "1k.DP" "norm" 4 &
#                wait_all
#              done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt
#              wait
#              while read line;do
#                sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/script_1M.sh $line
#              done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples.txt
#         fi 
#    fi
#done
#wait
#for num in {1..7}; do
#    for i in {A,B,D}; do
#        #1
#        for sample in {${SAMPLE1},${SAMPLE2}}; do
#            if [[ ! -s "tmp/${SAMPLE1}/chr${num}${i}_snp_${item}" ]];then
#                 if [[ -d "CNV_heatmap/${sample}" ]]&&[[ -s "CNV_heatmap/${sample}/chr${num}${i}.1M.norm" ]]; then
#                      ./muti_func_snp_compare.py -b 1000000 -i CNV_heatmap/${sample}/chr${num}${i}.1M.norm --mask_cnv on -o tmp/${sample}/chr${num}${i}.mask_CNV &
#                 else
#                      ./muti_func_snp_compare.py -b 1000000 -i /data2/rawdata2/readDepth/${sample}/chr${num}${i}.1M.norm --mask_cnv on -o tmp/${sample}/chr${num}${i}.mask_CNV &
#                 fi
#            fi
#        done
#        wait
#        #2
#        awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE2"_CNV"}' SAMPLE2="$SAMPLE2" tmp/${SAMPLE1}/chr${num}${i}.mask_CNV tmp/${SAMPLE2}/chr${num}${i}.mask_CNV > tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE2}to${SAMPLE1}_own_CNV
#        awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE1"_CNV"}' SAMPLE1="$SAMPLE1" tmp/${SAMPLE2}/chr${num}${i}.mask_CNV tmp/${SAMPLE1}/chr${num}${i}.mask_CNV > tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV
#        awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==1) print $1"\t"$2"\t"$3"\t""both_CNV"}' tmp/${SAMPLE1}/chr${num}${i}.mask_CNV tmp/${SAMPLE2}/chr${num}${i}.mask_CNV > tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}both_CNV
#        #3
#        cat tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE2}to${SAMPLE1}_own_CNV tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}both_CNV | sort -u -nk2,2 > tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV
#        #paste -d "\t" ${sample1}/chr${num}${i}.1M.norm ${sample2}/chr${num}${i}.1M.norm |sed '1d'|awk '{print $0"\t""chr""'$num'""'$i'"}' > chr${num}${i}_combine.1M.norm
#        #paste -d "\t" chr${num}${i}_combine.1M.norm ../lx99_jm22/chr${num}${i}.diff >> combine_gene_pos
#        wait_all
#    done
#done
#wait
#for num in {1..7};do
#    for i in {A,B,D};do
#        for item in {homo_DP,homo_GQ};do
#            if [[ ! -s "tmp/${SAMPLE1}/chr${num}${i}_snp_${item}" ]]; then    
#            ./muti_func_snp_filter.py -i ${total}/chr${num}${i}_snp.gt --single_DP_GQ on -s ${first} --chromsome tmp/${SAMPLE1}/chr${num}${i}.mask_CNV |grep ${item}|sort -nk1,1> tmp/${SAMPLE1}/chr${num}${i}_snp_${item} &
#            fi
#            if [[ ! -s "tmp/${SAMPLE2}/chr${num}${i}_snp_${item}" ]]; then
#            ./muti_func_snp_filter.py -i ${total}/chr${num}${i}_snp.gt --single_DP_GQ on -s ${second} --chromsome tmp/${SAMPLE2}/chr${num}${i}.mask_CNV |grep ${item}|sort -nk1,1> tmp/${SAMPLE2}/chr${num}${i}_snp_${item} &
#            fi
#            if [[ ! -s "tmp/${SAMPLE1}/chr${num}${i}_indel_${item}" ]]; then
#            ./muti_func_snp_filter.py -i ${total}/chr${num}${i}_indel.gt --single_DP_GQ on -s ${first} --chromsome tmp/${SAMPLE1}/chr${num}${i}.mask_CNV |grep ${item}|sort -nk1,1> tmp/${SAMPLE1}/chr${num}${i}_indel_${item} &   
#            fi
#            if [[ ! -s "tmp/${SAMPLE2}/chr${num}${i}_indel_${item}" ]]; then
#            ./muti_func_snp_filter.py -i ${total}/chr${num}${i}_indel.gt --single_DP_GQ on -s ${second} --chromsome tmp/${SAMPLE2}/chr${num}${i}.mask_CNV |grep ${item}|sort -nk1,1> tmp/${SAMPLE2}/chr${num}${i}_indel_${item} & 
#            fi
#            wait_all
#        done
#    done
#done
#wait
#
### level "1.5,2.85"
#for num in {1..7}; do
#    for i in {A,B,D};do
###   path="/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/"
#        ./muti_func_snp_compare.py -b 1000000 -i ${total}/chr${num}${i}_snp.gt --DP_low 3 --DP_high 20 --GQ_sample 8 -1 ${first} -2 ${second} --snp_graphy on --chromosome tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV -o tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}_${SAMPLE2}_unmatchhomo_snp_density &
#        ./muti_func_snp_compare.py --density_graph_raw_vcf_two on -i ${total}/chr${num}${i}_snp.gt -1 ${first} -2 ${second} --chromosome tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 --LEVEL "1.5,3.0" -o tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}_${SAMPLE2}unmatch_homo_snp_level &
##
#    done
#done
#wait
#cat tmp/${SAMPLE1}_${SAMPLE2}/chr*.${SAMPLE1}_${SAMPLE2}_unmatchhomo_snp_density |grep -v 'diff' > tmp/${SAMPLE1}_${SAMPLE2}/combine.${SAMPLE1}_${SAMPLE2}_unmatchhomo_snp_density
#./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${SAMPLE1}_${SAMPLE2}" -s ".${SAMPLE1}_${SAMPLE2}unmatch_homo_snp_level" > plotfile_${SAMPLE1}_${SAMPLE2}
#Rscript ~/R/graph.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${SAMPLE2} --sample2_name ${SAMPLE2_name}
for num in {1..7}; do
    for i in {A,B,D};do
        ./muti_func_snp_compare.py --transition_transversion_count on -i ${total}/chr${num}${i}_snp.gt --chromosome tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}_${SAMPLE2}unmatch_homo_snp_level  --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -1 ${first} -2 ${second} | sort -nk2,2 > tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_diff_transi_ver &
        ./muti_func_snp_compare.py --indel_count on -i ${total}/chr${num}${i}_indel.gt --chromosome tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}.${SAMPLE1}_${SAMPLE2}unmatch_homo_snp_level --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -1 ${first} -2 ${second}  | sort -nk2,2 > tmp/${SAMPLE1}_${SAMPLE2}/chr${num}${i}_indel_len_count &
    done
done
wait
cat tmp/${SAMPLE1}_${SAMPLE2}/chr*_diff_transi_ver | awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}/combine_total_diff_transi_ver
cat tmp/${SAMPLE1}_${SAMPLE2}/chr*_diff_transi_ver |grep -E 'high|mid'| awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}/combine_diff_diff_transi_ver
cat tmp/${SAMPLE1}_${SAMPLE2}/chr*_indel_len_count | awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}/combine_total_indel_len_count 
cat tmp/${SAMPLE1}_${SAMPLE2}/chr*_indel_len_count  |grep -E 'high|mid'| awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}/combine_diff_indel_len_count 

Rscript ~/R/transi_tranver.R -s ${SAMPLE1}_${SAMPLE2}
Rscript ~/R/indel_length.R -s ${SAMPLE1}_${SAMPLE2}
