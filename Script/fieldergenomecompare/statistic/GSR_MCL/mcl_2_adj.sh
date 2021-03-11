
mcl_file=$1
arritemidx(){
  local tmp 
  local count=0
  local array=`echo $1`
  for tmp in ${array[@]};do
    if test $2 = $tmp;then
      echo $count
      return
    fi  
    count=$(( $count + 1 ))
  done
  echo -1
}

arrslice(){
  array=($1)
  if [ $2 == -1 ];then
    echo ${array[@]}
  elif [ $2 == 0 ];then
    echo ${array[@]:1}
  else
    #echo "${array[@]:0:$2} ${array[@]:$(( $2 + 1 ))}"
    echo "${array[*]:0:$2} ${array[*]:$(( $2 + 1 ))}"
  fi
}
num_tmp=0
line_num=1


while read lines;do
arr=($lines)
arr_2=($lines)
for SAMPLE1 in ${arr[@]};do
    for SAMPLE2 in ${arr_2[@]};do
         sample1=`echo ${SAMPLE1} | awk -F '_' '{print $2}'`
         sample2=`echo ${SAMPLE2} | awk -F '_' '{print $2}'`
         num=`arritemidx "${arr_2[*]}" ${SAMPLE1}`
         arr_2=(`arrslice "${arr_2[*]}" $num`)
         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
             echo -e "${SAMPLE1}\t${SAMPLE2}\t1"
         fi
    done
#wait_all
line_num=$(( $line_num + 1 ))
done
done< $mcl_file
