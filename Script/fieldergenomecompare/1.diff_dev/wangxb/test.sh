#!/bin/bash
#while read line;do
#echo ${line} |awk '{split($0,a); len=length(a);for(i=1; i<=len; i++){split(a[i],b,"/");if(b[1]=="."){a[i]=6} else if (b[1]==b[2]&&b[1]==0){a[i]=0} else if (b[1]==b[2]){a[i]=1} else if (b[1]!=b[2]){a[i]=7}b[1];print a[i]}}'|tr '\n' '	'| sed s/'	'$/'\n'/g >> test.txt
#done < gsub_test.txt


#declare -A arr_2
#while read line;do
#arr_2[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $1}'`
#done < ${1}

#echo ${!arr_2[@]}
#Array.prototype.arrDel = function(n){
#    if (n<0){
#        return this;
#    }else{
#        return this.slice(0,n).concat(this.slice(n+1,this.length));
#    }
#}
BED=$1
for CHR in `cut -f1 $BED|sort | uniq`; do
n=`cat ${BED} | grep ${CHR} | wc -l`
echo $n $CHR
done
