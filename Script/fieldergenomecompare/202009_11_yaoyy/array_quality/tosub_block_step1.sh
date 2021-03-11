file_name=$1
if true;then
declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $3}'`]=`echo $line | awk '{print $2}'`
done < ../../metadata_cultivar_196_nongke_10+.txt

while read line;do
chr=`echo $line | awk '{print $1}'`
start=`echo $line | awk '{print $2}'`
sample=`echo $line | awk -F '[ _,]' '{print $4}'`
a=`bcftools view -v snps  -r ${chr}:${start} -s ${arr_3[${sample}]} ../../../08.mergeGVCF/196_9_8_dev/${chr}.bcf.gz|bcftools query -f "%REF\t%ALT\t%AF\t[%GT]\n"|awk -F '[\t/]' '{if ($4==0){print $1"\t"$2"\t"$3"\t"$1}else{print $1"\t"$2"\t"$3"\t"$2}}'`
echo $line $a | sed 's/ /\t/g' >> ${file_name}_maf_alt.txt
done <${file_name}.txt

while read line;do
chr=`echo $line | awk '{print $1}'`
start=`echo $line | awk '{print $2}'`
end=`echo $line | awk '{print $3}'`
block=`echo $line | awk '{print $5}'`
awk -vstart=$start -vend=$end -vblock=$block -vchr=$chr '{if ($2>=start && $2<=end && $1==chr) print $0"\t"block}' ${file_name}_maf_alt.txt >> ${file_name}_maf_alt_block.txt
done<gene_list_all_merge.txt
fi

n=`cut -f11 ${file_name}_maf_alt_block.txt| sort -n | tail -1`
num=0
for i in $(eval echo {1..${n}});do
awk -vi=${i} '{if ($11==i) print $0}' ${file_name}_maf_alt_block.txt | sort -k3,3 >tmp.txt

cut -f3  tmp.txt| uniq | awk -vnum=$num '{print $1"\t"NR+num}' > test.txt
line=`cat test.txt|wc -l `
num=$(( $num + $line ))
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next} NR!=FNR{if($3 in a ==1) print $11,$1,$2,$3,a[$3],$4,$5,$6,$7,$8,$9,$10}' test.txt tmp.txt >> ${file_name}_final.txt


done

