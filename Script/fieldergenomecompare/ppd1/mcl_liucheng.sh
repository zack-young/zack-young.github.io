#!/bin/bash
set -euxo pipefail
if false;then
for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do #chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#CHR="chr1A"
n=$(wc -l < ~/mapping/fieldergenomecompare/statistic/GSR_MCL/${CHR}.1M.bed)
parallel --xapply -j procfile ./bisample_compare_bin_edge_list.sh sample_list.txt ::: $(eval echo {1..$n}) ::: ${CHR} ::: $(eval echo {1..$n}) ::: homo_undefined_snp_level ::: /data/user/yangzz/mapping/fieldergenomecompare/10+_genome ::: 10+_GSR ::: ../metadata_cultivar_tmp.txt

done
wait
fi

if false;then
for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
n=$(wc -l < ~/mapping/fieldergenomecompare/statistic/GSR_MCL/${CHR}.1M.bed)
echo $n
(for i in `eval echo {1..$n}`;do
#grep -v '20' CNC_sample/${CHR}_homo_undefined_snp_level_${i} | sponge CNC_sample/${CHR}_homo_undefined_snp_level_${i}
/home/wangzh/bin/mcl 10+_GSR/${CHR}_homo_undefined_snp_level_${i} --abc -o 10+_mcl/${CHR}.${i}_mcl_undefined
#gsr=`cat mcl_dir/${CHR}.${i}_mcl_undefined | tr '\t' '\n' | wc -l | cut -d ' ' -f1`
cnv_del=`sed -n "${i}p"  /data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/CNV_count/${CHR}.combine_mask_CNV_deletion_count_compensent | cut -f3`
cnv_dup=`sed -n "${i}p"  /data/user/yangzz/mapping/fieldergenomecompare/CNV_stastistics/CNV_count/${CHR}.combine_mask_CNV_duplication_count_compensent | cut  -f3`
#if [ $clust == 0 ];then
#ratio=0
#else
#ratio=$(( $gsr/$clust ))
#fi
#./entropy.py -i CNL_mcl_dir/${CHR}.${i}_mcl_undefined -p 60 -c ${cnv_del} -d ${cnv_dup} >> ${CHR}_CNL_entropy
#./PIC.py -i CNL_mcl_dir/${CHR}.${i}_mcl_undefined -p 60 -c ${cnv_del} -d ${cnv_dup} >> ${CHR}_CNL_PIC
#./hap_number_count.py -i mcl_dir/${CHR}.${i}_mcl_undefined -p 198 -c ${cnv_del} -d ${cnv_dup} >> ${CHR}_all_hap
#num=$(( 198-$gsr+$clust-$cnv ))
#cat mcl_dir/${CHR}.${i}_mcl_undefined | tr '\t' '\n'  |awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($3 in a ==0) print $3}' - ../../metadata_cultivar_198.txt >> sample_stat
#for n in `eval echo {1..$clust}`;do
#num=`head -1 mcl_dir/${CHR}.${i}_mcl_undefined | tr '\t' '\n' | wc -l | cut -d ' ' -f1`
#echo -e "${CHR}\t$(( ($i-1)*1000000+1 ))\t$(( $i*1000000+1 ))\t${ratio}" >> ${CHR}_mcl_ratio_undefined_summary
#echo -e "${CHR}\t$(( ($i-1)*1000000+1 ))\t$(( $i*1000000+1 ))\t${num}" >> ${CHR}_mcl_num_undefined_summary_cnv
#echo -e "${CHR}\t$(( ($i-1)*1000000+1 ))\t$(( $i*1000000+1 ))\t${num}" >> ${CHR}_mcl_1st_clust_num_undefined_summary
done
) &
done
fi

if true;then
#for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do #chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
CHR="chr6A"
n=$(wc -l < ~/mapping/fieldergenomecompare/statistic/GSR_MCL/${CHR}.1M.bed)
parallel --xapply -j procfile ./bisample_compare_hap_type.sh ::: ${CHR} ::: $(eval echo {1..$n})

#done
wait
fi

