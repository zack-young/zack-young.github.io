#!/bin/bash
#set -euxo pipefail
CHR=$1
n=$2
typ=$3
#dev_path=$6
#out_dir=$7

START=$(( ($n-1)*1000000+1 ))

#$declare -A arr_3
#while read line;do
#arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
#done < ${meta_data}
#SAMPLE_list=${SAMPLE1}
#arr=(`echo ${SAMPLE_list} | tr ',' ' '`)
#arr_3=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
#declare -A arr_2
#while read line;do
#arr_2[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $1}'`
#done < ${meta_data}

num_tmp=0
line_num=1
if false;then
#grep -v '20' CN_sample/${CHR}_homo_undefined_snp_level_${i} | sponge CN_sample/${CHR}_homo_undefined_snp_level_${i}
awk '{if ($3!=20) {print $0}}' ${typ}_sample/${CHR}_homo_undefined_snp_level.sorted.v4_${n}  > ${typ}_sample/${CHR}_homo_undefined_snp_level.sorted.v4_noCNV_${n} 

/home/wangzh/bin/mcl ${typ}_sample/${CHR}_homo_undefined_snp_level.sorted.v4_noCNV_${n} --abc -o ${typ}_mcl_dir/${CHR}.${n}_mcl
#gsr=`cat mcl_dir/${CHR}.${i}_mcl_undefined | tr '\t' '\n' | wc -l | cut -d ' ' -f1`
fi
if false;then

dev='/data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics'
line_num=`awk -va=$START '{if($1==a) print $0}' ${dev}/CNV_${typ}_count/${CHR}.combine_mask_CNV_deletion_count_sample|wc -l `
if [ $line_num != 0 ];then
cnv_del=`awk -va=$START '{if($1==a) print $3}' ${dev}/CNV_${typ}_count/${CHR}.combine_mask_CNV_deletion_count_sample`
cnv_del_sample=`awk -va=$START '{if($1==a) print $4}' ${dev}/CNV_${typ}_count/${CHR}.combine_mask_CNV_deletion_count_sample`
else
cnv_del=0
cnv_del_sample='none'
fi

line_num=`awk -va=$START '{if($1==a) print $0}' ${dev}/CNV_${typ}_count/${CHR}.combine_mask_CNV_duplication_count_sample|wc -l`
if [ $line_num != 0 ];then
cnv_dup=`awk -va=$START '{if($1==a) print $3}' ${dev}/CNV_${typ}_count/${CHR}.combine_mask_CNV_duplication_count_sample`
cnv_dup_sample=`awk -va=$START '{if($1==a) print $4}' ${dev}/CNV_${typ}_count/${CHR}.combine_mask_CNV_duplication_count_sample`
else
cnv_dup=0
cnv_dup_sample='none'
fi
#echo $cnv_dup_sample $cnv_del_sample
./hap_number_count.py -i ALL_mcl_dir/${CHR}.${n}_mcl -r ${CHR} -s ${START} -c ${cnv_del}  -d ${cnv_dup} -t ${typ} --del_sample ${cnv_del_sample} --dup_sample ${cnv_dup_sample}
#echo ${cnv_dup} ${cnv_del}
num_C=`./entropy.py -i ALL_mcl_dir/${CHR}.${n}_mcl  -c ${cnv_del}  -d ${cnv_dup} -t ${typ}`
echo -e "${CHR}\t${START}\t${num_C}" >> ${CHR}_${typ}_entropy
fi

for i in CNL CNC;do
for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
for line in `cut -f2 ${chr}_ALL_hap_count.txt| sort -n | uniq` ; do
a=`grep $i ${chr}_ALL_hap_count.txt|grep -w $line | grep -vE 'exclusive|deletion|duplication' | wc -l`
b=`grep -w $line ${chr}_ALL_hap_count.txt| grep  exclusive|cut -f5 |tr ',' '\n'|grep  $i| wc -l`
echo $line	$(( $a+$b ))	$i
done
done
done
#ratio=0
#else
#ratio=$(( $gsr/$clust ))
#fi
#num_C=`./entropy.py -i CN_mcl_dir/${CHR}.${n}_mcl  -c ${cnv_del_C}  -d ${cnv_dup_C} -t C`
#echo ${CHR}	${START}	${num_C} >> ${CHR}_CNC_entropy
#num_L=`./entropy.py -i CN_mcl_dir/${CHR}.${n}_mcl  -c ${cnv_del_L}  -d ${cnv_dup_L} -t L`
#echo ${CHR}	${START}	${num_L} >> ${CHR}_CNL_entropy
#num_all=`./entropy.py -i CN_mcl_dir/${CHR}.${n}_mcl  -c ${cnv_del_all}  -d ${cnv_dup_all} -t all`
#echo ${CHR}	${START}	${num_all} >> ${CHR}_CN_entropy
#./PIC.py -i CNL_mcl_dir/${CHR}.${i}_mcl_undefined -p 60 -c ${cnv_del} -d ${cnv_dup} >> ${CHR}_CNL_PIC
#num=$(( 198-$gsr+$clust-$cnv ))
#cat mcl_dir/${CHR}.${i}_mcl_undefined | tr '\t' '\n'  |awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($3 in a ==0) print $3}' - ../../metadata_cultivar_198.txt >> sample_stat
#for n in `eval echo {1..$clust}`;do
#num=`head -1 mcl_dir/${CHR}.${i}_mcl_undefined | tr '\t' '\n' | wc -l | cut -d ' ' -f1`
#echo -e "${CHR}\t$(( ($i-1)*1000000+1 ))\t$(( $i*1000000+1 ))\t${ratio}" >> ${CHR}_mcl_ratio_undefined_summary
#echo -e "${CHR}\t$(( ($i-1)*1000000+1 ))\t$(( $i*1000000+1 ))\t${num}" >> ${CHR}_mcl_num_undefined_summary_cnv
#echo -e "${CHR}\t$(( ($i-1)*1000000+1 ))\t$(( $i*1000000+1 ))\t${num}" >> ${CHR}_mcl_1st_clust_num_undefined_summarydone
