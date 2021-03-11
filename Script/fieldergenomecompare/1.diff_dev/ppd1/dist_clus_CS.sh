#!/bin/bash

CHR=$1
N=$2
snp_path=$3
sample_list=$4
#bcftools query -f "[%GT\t]\n" -H -R <(tail -n +${N} ${CHR}.1M.bed|head -n 1) /data2/rawdata2/mergeFile/mergeTetra/all/${CHR}.ann.bcf.gz -S ../origin/data/WGS_samplelist.txt | sed 's$0/0$0$g;s$1/1$2$g;s$0/1$1$g;s$./.$6$g;s/\t$//' > ${CHR}.${N}.rawgeno
#python cluster.py -i ${CHR}.${N}.rawgeno -o ${CHR}.${N}
#rm ${CHR}.${N}.rawgeno

# tail -n +a  read file from a line to the end

#===========================yangzz version
bcftools view -v snps --min-ac=1 -M2 -m2 -R <(tail -n +${N} ${CHR}.1M.bed|head -n 1) /data/user/chenym/100rna-seq/guoyw_leaf_root_h2o2_heterogeneous_system/00.reseq/02.map/total.filter.final.merge.bcf.gz -S ${sample_list} > ${CHR}.${N}.rawgeno
bcftools query -f "[%GT\t]\n" ${CHR}.${N}.rawgeno >  ${CHR}.${N}.rawgt

touch ${CHR}.${N}.matrix
while read line;do
echo ${line} |awk '{split($0,a); len=length(a);for(i=1; i<=len; i++){split(a[i],b,"/");if(b[1]=="."){a[i]=30} else if (b[1]==b[2]&&b[1]==0){a[i]=1} else if (b[1]==b[2]){a[i]=5} else if (b[1]!=b[2]){a[i]=30};print a[i]}}'|tr '\n' '	'| sed s/'	'$/'\n'/g >>  ${CHR}.${N}.matrix
done < ${CHR}.${N}.rawgt

cat header.txt ${CHR}.${N}.matrix|sponge ${CHR}.${N}.matrix

bcftools query -f "[%DP\t]\n" -H ${CHR}.${N}.rawgeno |sed 's/\./0/g'> ${CHR}.${N}.DPmatrix


bcftools query -f "[%GQ\t]\n" -H ${CHR}.${N}.rawgeno |sed 's/\./0/g'> ${CHR}.${N}.GQmatrix

python cluster_yangzz_CS.py -i ${CHR}.${N}.matrix -d ${CHR}.${N}.DPmatrix -g ${CHR}.${N}.GQmatrix -u ${N} -n /data/user/yangzz/mapping/fieldergenomecompare/5.cross_sample/CNV_masker/${CHR}.CNVfilter.txt -o ${CHR}.${N}.diff

rm ${CHR}.${N}.rawgt
rm ${CHR}.${N}.rawgeno
rm ${CHR}.${N}.matrix
rm ${CHR}.${N}.DPmatrix
rm ${CHR}.${N}.GQmatrix


