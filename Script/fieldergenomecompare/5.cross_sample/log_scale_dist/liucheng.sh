array1=()
sample_compare_path="/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare"
declare -A dic
#追加字典
num=1
while read line;do
a=`echo $line | awk '{print$1}'`
dic+=([${a}]="${num}")
num=$(( $num + 1 ))
done < /data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_final.txt
if true; then
while read line;do
first=`echo $line |awk '{print $1}'`
for i in `echo $line |awk '{ for (i=2; i<=10; i++) print $i }'|tr '\n' '\t'|sed s/"\t"$//g`;do #change 3 to NF to print 10
if [[ -d "${sample_compare_path}/${first}_${i}" ]]; then
array1[${#array1[@]}]=${dic["${first}"]}'_'${dic["${i}"]}
else
array1[${#array1[@]}]=${dic["${i}"]}'_'${dic["${first}"]}
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

