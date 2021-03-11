#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#CHR='chr3A'
#while read line;do
#num=`echo ${line}| awk '{print $1"\t"$2"\t"$3}'  |bedtools intersect -b - -a gene_list_all.txt | wc -l`
#echo $line	$num
#if [ $num != 0 ];then
#echo ${line}| awk '{print $1"\t"$2"\t"$3}' | bedtools intersect -b - -a gene_list_all.txt |bedtools merge -i - -d 1000000 -c 4 -o distinct
#fi
#done < ../../1.diff_dev/${CHR}.1M.bed
#done
#awk '{print $1"\t"int($2/1000000)*1000000"\t"int($2/1000000)*1000000+1"\t"$4}' gene_list_all.txt | bedtools merge -i - -c 4 -o distinct > gene_list_all_merge.txt
#gawk -F '[,\t]' 'ARGIND==1{a[$1]=$1} ARGIND==2 {k=0;for (i=3;i<=NF;i++){if($i in a){k++}}print $0"\t"k}' good_q_sample.txt gene_list_hap_type.txt > gene_list_hap_type_good_count.txt
#awk '{if($6!=0)print $0}' gene_list_hap_type_good_count.txt|grep -Ev 'deletion|duplication'|sort -k1,1 -k2,2n > gene_list_hap_type_good_count_filter.txt

#while read line;do
#pos=`echo $line|awk '{print $1$2}'`
#echo $line|awk '{print $1"\t"$2"\t"$3}' >> gene_list_N_hap.txt
#echo $line|awk '{print $4}' >> gene_list_N_hap.txt
#awk '{if ($1$2==pos) {print $3"\t"$4"\t"$5"\t"$6}}' pos="$pos"  gene_list_hap_type_good_count.txt >> gene_list_N_hap.txt #combine_all_hap.txt
#done<gene_list_all_merge.txt


#awk  'BEGIN{OFS="\t"}{print $1,$2,$2,$5}' test.txt | bedtools merge -i     - -c 4 -o distinct -delim "|"
#./hap_number_count.py -i gene_list_hap_type_good_count_filter.txt|awk  'BEGIN{OFS="\t"}{print $1,$2,$2+1,$5}' | bedtools merge -i     - -c 4 -o distinct -delim "|" >gene_list_hap_type_good_count_filter_merge.txt
for i in  $(eval echo {1..2354673..50000});do
next=$(( $i + 49999))
parallel -j procfile sh tosub_inside.sh ::: chr1B ::: 1 ::: 239000000 ::: $(eval echo {$i..$next}) ::: ../../metadata_cultivar_196_nongke_10+.txt ::: 'CNC_NongDa3338,CNC_JingDong6,CNC_Shi4185,CNC_YuMai21,EUC_Lovrin10,CNC_ChangZhi6406,CNC_NongDa5181,CNC_ZhouMai18,CNC_LuMai15,CNC_JiNan16,CNC_ZhongMai875,CNC_HengGuan35,CNC_ZhouMai22,CNC_ZhouMai16,CNC_LunXuan987,CNC_Ak58,CNC_Aimengniu,CNC_Een1,CNC_Xumai856' ::: $i
done

