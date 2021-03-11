while read line;do
a=`echo $line | cut -f2`
b=`grep $a /data/annotation/wheat/CS_IWGSC/v1p1/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3 | grep gene|cut -f1,4,5`
echo -e "$line\t$b" >> gene_list_new.txt\ndone<new_gene_list.txt

awk '{print $1"\t"int($2/1000000)*1000000"\t"int($2/1000000)*1000000+1"\t"$4}' gene_list_all.txt | bedtools merge -i - -c 4 -o distinct > gene_list_all_merge.txt
gawk 'ARGIND==1{a[$1$2]=$0} ARGIND==2{if($1$2 in a){print $0}}' gene_list_all_merge.txt combine_all_hap.txt > gene_list_hap_type.txt
#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#while read line;do
#num=`echo ${line}| awk '{print $1"\t"$2"\t"$3}'  |bedtools intersect -b - -a gene_list_all.txt | wc -l`
#echo $line	$num
#if [ $num != 0 ];then
#echo ${line}| awk '{print $1"\t"$2"\t"$3}' | bedtools intersect -b - -a gene_list_all.txt |bedtools merge -i - -d 1000000 -c 4 -o distinct >> gene_list_N_hap.txt
#fi
#done < ../../1.diff_dev/${CHR}.1M.bed
#done
gawk -F '[,\t]' 'ARGIND==1{a[$1]=$1} ARGIND==2 {k=0;for (i=3;i<=NF;i++){if($i in a){k++}}print $0"\t"k}' good_q_sample.txt gene_list_hap_type.txt > gene_list_hap_type_good_count.txt
awk '{if($6!=0)print $0}' gene_list_hap_type_good_count.txt|grep -Ev 'deletion|duplication'|sort -k1,1 -k2,2n > gene_list_hap_type_good_count_filter.txt

#while read line;do
#pos=`echo $line|awk '{print $1$2}'`
#echo $line|awk '{print $1"\t"$2"\t"$3}' >> gene_list_N_hap.txt
#echo $line|awk '{print $4}' >> gene_list_N_hap.txt
#awk '{if ($1$2==pos) {print $3"\t"$4"\t"$5"\t"$6}}' pos="$pos"  gene_list_hap_type_good_count.txt >> gene_list_N_hap.txt #combine_all_hap.txt
#done<gene_list_all_merge.txt


awk  'BEGIN{OFS="\t"}{print $1,$2,$2,$5}' test.txt | bedtools merge -i     - -c 4 -o distinct -delim "|"
./hap_number_count.py -i gene_list_hap_type_good_count_filter.txt|awk  'BEGIN{OFS="\t"}{print $1,$2,$2+1,$5}' | bedtools merge -i     - -c 4 -o distinct -delim "|" >gene_list_hap_type_good_count_filter_merge.txt

