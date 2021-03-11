#!/bin/bash
set -eu  
#uxo pipefail
#trying to find genome shared proportion of each accessions
#sample=$1
dev1="/data/user/yangzz/mapping/fieldergenomecompare"
dev="/data3/user3/wangwx/projs/HMM_for_yzz_comp/201212/201_sample_compare_HMMdata"
if false;then
for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do       
#sort  -nk 2,2 ${chr}.homo_hmm_snp_level |cut -f4|sed '1i\LEVEL'  |paste ${chr}.1M.density -|sed s/low/sim/g|sed s/high/diff/g|sed s/mid/diff/g > ${chr}.1M.density_level
#chr="chr1A"
(while read line;do
item=`echo $line |awk '{print $1}'`
name=`echo $line |awk '{print $3}'`
lis_CNC_CNL=`grep ${item} sample_path_CNC_CNL.txt| awk '{print dev"/"$1"/"chr".homo_undefined_snp_level.sorted.v4"}' dev="$dev" chr="$chr"| tr '\n' ' '|sed s/' '$//g`
lis_CNC_FC=`grep ${item} sample_path_CNC_FC.txt| awk '{print dev"/"$1"/"chr".homo_undefined_snp_level.sorted.v4"}' dev="$dev" chr="$chr"| tr '\n' ' '|sed s/' '$//g`
num_sim_L=`cat ${lis_CNC_CNL} | grep -E 'low'|cut -f1,2|sort -k1,2|uniq|wc  -l`
num_cnv_L=`cat ${lis_CNC_CNL} | grep -E 'both'|cut -f1,2|sort -k1,2|uniq|wc  -l`
num_sim_F=`cat ${lis_CNC_FC} | grep -E 'low'|cut -f1,2|sort -k1,2|uniq|wc  -l`
num_cnv_F=`cat ${lis_CNC_FC} | grep -E 'both'|cut -f1,2|sort -k1,2|uniq|wc  -l`

echo $num_sim_L $num_cnv_L $num_sim_F $num_cnv_F $name >> ${chr}.CNC2CNL_FC_sim_count.txt

done < sample_list_CNC.txt ) &
done
wait
fi
if true;then
while read line;do
name=`echo $line | awk '{print $3}'`
num=`cat chr*.CNC2CNL_FC_sim_count.txt | grep -w ${name} | awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'`
num2=`cat chr*.CNC2CNL_FC_sim_count.txt | grep -w ${name} | awk 'BEGIN{sum=0}{sum+=$2}END{print sum}'`
num3=`cat chr*.CNC2CNL_FC_sim_count.txt | grep -w ${name} | awk 'BEGIN{sum=0}{sum+=$3}END{print sum}'`
num4=`cat chr*.CNC2CNL_FC_sim_count.txt | grep -w ${name} | awk 'BEGIN{sum=0}{sum+=$4}END{print sum}'`
echo  $num $num2 $num3 $num4 $name >> CNC_sample_CNL_FC_sim_count.txt
done < sample_list_CNC.txt
fi
