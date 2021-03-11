#!/bin/bash

CHR="chr1B"
N="S"
#snp_path=$3
#sample_list=$4
#bcftools query -f "[%GT\t]\n" -H -R <(tail -n +${N} ${CHR}.1M.bed|head -n 1) /data2/rawdata2/mergeFile/mergeTetra/all/${CHR}.ann.bcf.gz -S ../origin/data/WGS_samplelist.txt | sed 's$0/0$0$g;s$1/1$2$g;s$0/1$1$g;s$./.$6$g;s/\t$//' > ${CHR}.${N}.rawgeno
#python cluster.py -i ${CHR}.${N}.rawgeno -o ${CHR}.${N}
#rm ${CHR}.${N}.rawgeno

# tail -n +a  read file from a line to the end

#===========================yangzz version
bcftools view -v snps --min-ac=1 -M2 -m2 -r chr1B:1-239000000   /data2/rawdata2/shinyKiloProj/KiloProj.ann.bcf.gz  > ${CHR}.${N}.rawgeno
bcftools query -f "[%GT\t]\n" ${CHR}.${N}.rawgeno >  ${CHR}.${N}.rawgt

touch ${CHR}.${N}.matrix
while read line;do
echo ${line} |awk '{split($0,a); len=length(a);for(i=1; i<=len; i++){split(a[i],b,"/");if(b[1]=="."){a[i]=30} else if (b[1]==b[2]&&b[1]==0){a[i]=1} else if (b[1]==b[2]){a[i]=5} else if (b[1]!=b[2]){a[i]=2};print a[i]}}'|tr '\n' '	'| sed s/'	'$/'\n'/g >>  ${CHR}.${N}.matrix
done < ${CHR}.${N}.rawgt

#cat header.txt ${CHR}.${N}.matrix|sponge ${CHR}.${N}.matrix

#python cluster_yangzz.py -i ${CHR}.${N}.matrix -d ${CHR}.${N}.DPmatrix -g ${CHR}.${N}.GQmatrix -u ${N} -n /data/user/yangzz/mapping/fieldergenomecompare/20200418_shannong/${CHR}.CNVfilter.txt -o ${CHR}.${N}.diff


