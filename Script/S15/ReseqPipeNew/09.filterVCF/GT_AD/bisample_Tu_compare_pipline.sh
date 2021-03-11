#!/bin/bash
set -x

while getopts 'c:v:a:s:m:b:d:f:n:' opt
do
    case "$opt" in
    c) SAMPLE1=$OPTARG ;;
    v) SAMPLE1_name=$OPTARG ;;
    a) SAMPLE1_vcf=$OPTARG ;;
    s) SAMPLE2=$OPTARG ;;
    m) SAMPLE2_name=$OPTARG ;;
    b) SAMPLE2_vcf=$OPTARG ;;
    d) total=$OPTARG ;;
    f) first=$OPTARG ;;
    n) second=$OPTARG ;;
    esac
done
time_date="`date +%Y-%m-%d%H%M%S`" 
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
if [[ -z "${SAMPLE1_vcf}" ]]; then
    echo "No second sample vcf name provided."; exit 1
fi
if [[ -z "${SAMPLE2_vcf}" ]]; then
    echo "No second sample vcf name provided."; exit 1
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
if [[ ! -d "tmp/${SAMPLE1}_${time_date}" ]]; then
     mkdir tmp/${SAMPLE1}_${time_date}
fi
if [[ ! -d "tmp/${SAMPLE2}_${time_date}" ]]; then
     mkdir tmp/${SAMPLE2}_${time_date}
fi
if [[ ! -d "tmp/${SAMPLE1}_${SAMPLE2}_${time_date}" ]]; then
     mkdir tmp/${SAMPLE1}_${SAMPLE2}_${time_date}
fi
if [[ ! -d "${SAMPLE1}_${SAMPLE2}" ]]; then
     mkdir ${SAMPLE1}_${SAMPLE2}
fi
#for VCF in /data/user/yangzz/mapping/S15/ReseqPipeNew/08.mergeGVCF/${total}/Tu?.bcf.gz; do
#for VCF in /data2/rawdata2/mergeFile/version2/chr??.bcf.gz; do
#for VCF in /data/user/yangzz/mapping/08.mergeGVCF/field_cultivar/chr??.ann.bcf.gz; do # only for cultivar and fielder
    #CHR=`basename $VCF | sed s/.ann.bcf.gz//g`
#    CHR=`basename $VCF | sed s/.bcf.gz//g`
#    echo "${VCF}"
#    if [[ ! -s "${SAMPLE1}_${SAMPLE2}/${CHR}_snp.gt" ]]; then
#         bcftools view  -v snps -s ${SAMPLE1_vcf},${SAMPLE2_vcf} -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > ${SAMPLE1}_${SAMPLE2}/${CHR}_snp.gt &
#    fi
#    if [[ ! -s "${SAMPLE1}_${SAMPLE2}/${CHR}_indel.gt" ]]; then
#         bcftools view  -v indels -s ${SAMPLE1_vcf},${SAMPLE2_vcf} -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > ${SAMPLE1}_${SAMPLE2}/${CHR}_indel.gt &
#    fi
#    wait_all
#done
#wait
for a in {${SAMPLE1},${SAMPLE2}}; do
    if [[ ! -d "field_cultivar/${a}" ]]&&[[ ! -d "CNV_heatmap/${a}_${time_date}" ]]; then 
             
              mkdir CNV_heatmap/${a}_${time_date}
              path_SM="/data2/rawdata2/readDepth_ori_1k/${a}"
              # norm       
              for chr in Tu1 Tu2 Tu3 Tu4 Tu5 Tu6 Tu7;do
              gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' $path_SM/${chr}.1.1k.DP > CNV_heatmap/${a}_${time_date}/${chr}.1.1M.tmp
              gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' $path_SM/${chr}.2.1k.DP > CNV_heatmap/${a}_${time_date}/${chr}.2.1M.tmp
              done
              cat CNV_heatmap/${a}_${time_date}/Tu*.1M.tmp > CNV_heatmap/${a}_${time_date}/combine_1M_DP
              sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/script_Tu_1M.sh CNV_heatmap/${a}_${time_date} &
    fi 
done
wait
for chr in Tu1 Tu2 Tu3 Tu4 Tu5 Tu6 Tu7; do
        #1
        for sample in {${SAMPLE1},${SAMPLE2}}; do
                 #if [[ -d "CNV_heatmap/${sample}" ]]&&[[ -s "CNV_heatmap/${sample}/chr${num}${i}.1M.norm" ]]; then
                 #     ./muti_func_snp_compare.py -b 1000000 -i CNV_heatmap/${sample}/chr${num}${i}.1M.norm --mask_cnv on -o tmp/${sample}_${time_date}/chr${num}${i}.mask_CNV &
                 if [[ -d "field_cultivar/${sample}" ]]&&[[ -s "field_cultivar/${sample}/${chr}.1M.norm" ]]; then
                      ./muti_func_snp_compare.py -b 1000000 -i field_cultivar/${sample}/${chr}.1M.norm --mask_cnv on -o tmp/${sample}_${time_date}/${chr}.mask_CNV 
                 else
                      ./muti_func_snp_compare.py -b 1000000 -i CNV_heatmap/${sample}_${time_date}/${chr}.1M.norm --mask_cnv on -o tmp/${sample}_${time_date}/${chr}.mask_CNV
                      #./muti_func_snp_compare.py -b 1000000 -i /data2/rawdata2/readDepth/${sample}/chr${num}${i}.1M.norm --mask_cnv on -o tmp/${sample}_${time_date}/chr${num}${i}.mask_CNV &
                 fi
        done
        #2
        for item in {deletion,duplication};do
            for sample in {${SAMPLE1},${SAMPLE2}}; do
                grep "${item}" tmp/${sample}_${time_date}/${chr}.mask_CNV > tmp/${sample}_${time_date}/${chr}.mask_CNV_${item}
            done
            a=`du -h tmp/${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item}`
            b=`echo $a |awk '{print $1}'`
            if [ "${b}" == 0 ];then
                echo 'file zero'
                awk '{print $1"\t"$2"\t"$3"\t"SAMPLE2"'"$item"'""_CNV"}' SAMPLE2="$SAMPLE2" tmp/${SAMPLE2}_${time_date}/${chr}.mask_CNV_${item}  > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item}
            else
                echo 'file not zero'
                awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE2"'"$item"'""_CNV"}' SAMPLE2="$SAMPLE2" tmp/${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item} tmp/${SAMPLE2}_${time_date}/${chr}.mask_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item}
            fi
            #
            a=`du -h tmp/${SAMPLE2}_${time_date}/${chr}.mask_CNV_${item}`
            b=`echo $a |awk '{print $1}'`
            echo ${b}
            if [ "${b}" == 0 ];then
                echo 'file zero'
                awk '{print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" tmp/${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item}  > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
            else
                echo 'file not zero'
                awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" tmp/${SAMPLE2}_${time_date}/${chr}.mask_CNV_${item} tmp/${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
            fi
            awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==1) print $1"\t"$2"\t"$3"\t""'"$item"'""_both_CNV"}' tmp/${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item} tmp/${SAMPLE2}_${time_date}/${chr}.mask_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}_${SAMPLE1}to${SAMPLE2}both_CNV_${item}
        #3
            cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item} tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item} tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}_${SAMPLE1}to${SAMPLE2}both_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV_${item}
        #paste -d "\t" ${sample1}/chr${num}${i}.1M.norm ${sample2}/chr${num}${i}.1M.norm |sed '1d'|awk '{print $0"\t""chr""'$num'""'$i'"}' > chr${num}${i}_combine.1M.norm
        #paste -d "\t" chr${num}${i}_combine.1M.norm ../lx99_jm22/chr${num}${i}.diff >> combine_gene_pos
done
        #mv tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV_deletion  tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}.${SAMPLE1}to${SAMPLE2}all_CNV
        cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV_deletion tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV_duplication | sort -u -nk2,2 > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV
    done
#
#for chr in Tu1 Tu2 Tu3 Tu4 Tu5 Tu6 Tu7; do
#        for item in {homo_DP,homo_GQ};do
#            if [[ ! -s "tmp/${SAMPLE1}_${time_date}/${chr}_snp_${item}" ]]; then    
#            ./muti_func_snp_filter.py -i ${SAMPLE1}_${SAMPLE2}/${chr}_snp.gt --single_DP_GQ on -s ${first} --chromsome tmp/${SAMPLE1}_${time_date}/${chr}.mask_CNV |grep ${item}|sort -nk1,1> tmp/${SAMPLE1}_${time_date}/${chr}_snp_${item} &
#            fi
#            if [[ ! -s "tmp/${SAMPLE2}_${time_date}/${chr}_snp_${item}" ]]; then
#            ./muti_func_snp_filter.py -i ${SAMPLE1}_${SAMPLE2}/${chr}_snp.gt --single_DP_GQ on -s ${second} --chromsome tmp/${SAMPLE2}_${time_date}/${chr}.mask_CNV |grep ${item}|sort -nk1,1> tmp/${SAMPLE2}_${time_date}/${chr}_snp_${item} &
#            fi
#            if [[ ! -s "tmp/${SAMPLE1}_${time_date}/${chr}_indel_${item}" ]]; then
#            ./muti_func_snp_filter.py -i ${SAMPLE1}_${SAMPLE2}/${chr}_indel.gt --single_DP_GQ on -s ${first} --chromsome tmp/${SAMPLE1}_${time_date}/${chr}.mask_CNV |grep ${item}|sort -nk1,1> tmp/${SAMPLE1}_${time_date}/${chr}_indel_${item} &   
#            fi
#            if [[ ! -s "tmp/${SAMPLE2}_${time_date}/${chr}_indel_${item}" ]]; then
#            ./muti_func_snp_filter.py -i ${SAMPLE1}_${SAMPLE2}/${chr}_indel.gt --single_DP_GQ on -s ${second} --chromsome tmp/${SAMPLE2}_${time_date}/${chr}.mask_CNV |grep ${item}|sort -nk1,1> tmp/${SAMPLE2}_${time_date}/${chr}_indel_${item} & 
#            fi
#            wait_all
#        done
#done
#wait

## level "1.5,2.85"
for chr in Tu1 Tu2 Tu3 Tu4 Tu5 Tu6 Tu7; do
        #./tu_muti_func_snp_compare.py -b 1000000 -i ${SAMPLE1}_${SAMPLE2}/${chr}_snp.gt --DP_low 3 --DP_high 99 --GQ_sample 8 -1 ${first} -2 ${second} --snp_graphy on --chromosome tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV -o tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}.${SAMPLE1}_${SAMPLE2}_unmatchhomo_snp_density &
        ./tu_muti_func_snp_compare.py --density_graph_raw_vcf_two on -i ${SAMPLE1}_${SAMPLE2}/${chr}_snp.gt -1 ${first} -2 ${second} --chromosome tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV --DP_low 3 --DP_high 99 --GQ_sample 8 -b 100000 --LEVEL "1.5,3.0" -o tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/${chr}.${SAMPLE1}_${SAMPLE2}unmatch_homo_snp_level &
#
done
wait
cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/Tu*.${SAMPLE1}_${SAMPLE2}_unmatchhomo_snp_density |grep -v 'diff' > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/combine.${SAMPLE1}_${SAMPLE2}_unmatchhomo_snp_density
./ptf_maker_Tu.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${SAMPLE1}_${SAMPLE2}_${time_date}" -s ".${SAMPLE1}_${SAMPLE2}unmatch_homo_snp_level" > plotfile/plotfile_${SAMPLE1}_${SAMPLE2}
./ptf_maker_Tu.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${SAMPLE1}_${SAMPLE2}_${time_date}" -s ".${SAMPLE1}to${SAMPLE2}all_CNV_deletion" > plotfile/plotfile_${SAMPLE1}_${SAMPLE2}_deletion
./ptf_maker_Tu.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${SAMPLE1}_${SAMPLE2}_${time_date}" -s ".${SAMPLE1}to${SAMPLE2}all_CNV_duplication" > plotfile/plotfile_${SAMPLE1}_${SAMPLE2}_duplication
Rscript ~/R/graph_diff_sim_Tu.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${SAMPLE2} --sample2_name ${SAMPLE2_name} -s ${time_date}
#Rscript ~/R/graph_CNV_deletion.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${SAMPLE2} --sample2_name ${SAMPLE2_name} -s ${time_date}
#Rscript ~/R/graph_CNV_duplication.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${SAMPLE2} --sample2_name ${SAMPLE2_name} -s ${time_date}

#Rscript ~/R/graph_CNV.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${SAMPLE2} --sample2_name ${SAMPLE2_name} -s ${time_date}
#for num in {1..7}; do
#    for i in {A,B,D};do
#        ./muti_func_snp_compare.py --transition_transversion_count on -i ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_snp.gt --chromosome tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}.${SAMPLE1}_${SAMPLE2}unmatch_homo_snp_level  --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -1 ${first} -2 ${second} | sort -nk2,2 > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_diff_transi_ver &
#        ./muti_func_snp_compare.py --indel_count on -i ${SAMPLE1}_${SAMPLE2}/chr${num}${i}_indel.gt --chromosome tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}.${SAMPLE1}_${SAMPLE2}unmatch_homo_snp_level --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 -1 ${first} -2 ${second}  | sort -nk2,2 > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_indel_len_count &
#    done
#done
#wait
#cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr*_diff_transi_ver | awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/combine_total_diff_transi_ver
#cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr*_diff_transi_ver |grep -E 'high|mid'| awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/combine_diff_diff_transi_ver
#cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr*_indel_len_count | awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/combine_total_indel_len_count 
#cat tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr*_indel_len_count  |grep -E 'high|mid'| awk  '{arr[$1] = arr[$1] + $2}END{for (a in arr) print a"\t"arr[a]}' > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/combine_diff_indel_len_count 
#
#Rscript ~/R/transi_tranver.R -s ${SAMPLE1}_${SAMPLE2}_${time_date}
#Rscript ~/R/indel_length.R -s ${SAMPLE1}_${SAMPLE2}_${time_date}
