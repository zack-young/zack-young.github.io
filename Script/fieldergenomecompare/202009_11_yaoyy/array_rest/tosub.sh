mergen_sample_lis=$1
ord=$2

meta="/data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_196_nongke_10+.txt"
sample_list=`cut -f2 ${meta}|tr '\n' ','|sed  s/','$/''/g`
vcf_path="/data/user/yangzz/mapping/08.mergeGVCF/196_9_8_dev"

bcftools query -f "[%GT\t]\n"  -H  ${vcf_path}/chr7D.bcf.gz -s  ${sample_list}|head -1 > header.txt
while read line;do
sample_list=`cut -f2 ${meta}|tr '\n' ','|sed  s/','$/''/g`
CHR=`echo ${line} | awk '{print $1}'`
START=`echo ${line} | awk '{print $2}'`
END=`echo ${line} | awk '{print $2+1000000}'`
com_sample=`echo ${line} | awk '{print $4}'`
##---------------
bcftools view -v snps -q 0.01:minor -s ${sample_list} -r ${CHR}:${START}-${END} --min-ac=1 -M2 -m2  ${vcf_path}/${CHR}.bcf.gz  > rawdata/${CHR}.${START}-${END}.rawgeno
bcftools query -f "[%GT\t]\n" -s ${sample_list} rawdata/${CHR}.${START}-${END}.rawgeno > rawdata/${CHR}.${START}-${END}.rawgt
bcftools query -f "%POS\n" rawdata/${CHR}.${START}-${END}.rawgeno >  rawdata/${CHR}.${START}-${END}.POS
line=`wc -l rawdata/${CHR}.${START}-${END}.rawgt| cut -d ' ' -f1`
if [[ "$line" != 0 ]];then
touch rawdata/${CHR}.${START}-${END}.matrix
while read line;do
echo ${line} |awk '{split($0,a); len=length(a);for(i=1; i<=len; i++){split(a[i],b,"/");if(b[1]=="."){a[i]=30} else if (b[1]==b[2]&&b[1]==0){a[i]=1} else if (b[1]==b[2]){a[i]=5} else if (b[1]!=b[2]){a[i]=11};print a[i]}}'|tr '\n' '\t'| sed s/'\t'$/'\n'/g >>  rawdata/${CHR}.${START}-${END}.matrix
done < rawdata/${CHR}.${START}-${END}.rawgt

cat header.txt rawdata/${CHR}.${START}-${END}.matrix|sponge rawdata/${CHR}.${START}-${END}.matrix
bcftools query -f "[%DP\t]\n" -s ${sample_list} -H rawdata/${CHR}.${START}-${END}.rawgeno |sed 's/\./0/g'> rawdata/${CHR}.${START}-${END}.DPmatrix


bcftools query -f "[%GQ\t]\n" -s ${sample_list} -H rawdata/${CHR}.${START}-${END}.rawgeno |sed 's/\./0/g'> rawdata/${CHR}.${START}-${END}.GQmatrix
line_num=`wc -l rawdata/${CHR}.${START}-${END}.POS| cut -d ' ' -f1`
parallel -j procfile sh tosub_inside.sh ::: ${CHR} ::: ${START} ::: ${END} ::: $(eval echo {1..$line_num}) ::: ${meta} ::: ${com_sample} ::: ${ord}
fi
##--------
done < ${mergen_sample_lis}
