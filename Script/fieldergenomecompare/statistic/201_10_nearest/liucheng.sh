array1=()
sample_compare_path="/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare"
num=1
if true; then
while read line;do
first=`echo $line |awk '{print $1}'`
for i in `echo $line |awk '{ for (i=2; i<=3; i++) print $i }'|tr '\n' '\t'|sed s/"\t"$//g`;do #change 3 to NF to print 10
if [[ -d "${sample_compare_path}/${first}_${i}" ]]; then
array1[${#array1[@]}]=${first}_${i}
else
array1[${#array1[@]}]=${i}_${first}
num+=1
fi
done
done < nearest10.txt
#echo ${array1[*]} 
#echo ${#array1[@]}
arr=($(echo ${array1[*]}|sed 's/ /\n/g'|sort | uniq))
#echo ${arr[*]}|awk '{ for (i=1; i<=NF; i++) print "/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare/"$i }'
echo ${arr[*]}|awk '{ for (i=1; i<=NF; i++) print $i }'
fi




if false;then
for chr in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
while read line;do
echo $line | awk '{ for (i=4; i<=NF; i += 4) print $i }' | tr '\n' '\t'|sed s/"\t"$/'\n'/g >> ${chr}.2_col_density
done < ${chr}.2_col_combine_density
done
fi
if false;then
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
cat `awk '{print $1"/"CHR".1M.delcnv_hete_density"}' CHR="$CHR" sample_path.txt | tr '\n' ' '`> ${CHR}.1_hete_combine_density &
wait_all
done
fi
