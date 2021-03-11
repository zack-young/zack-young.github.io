#funWithReturn(){
#  echo "这个函数会对输入的两个数字进行相加运算..."
#  echo "输入第一个数字: "
#  read aNum
#  echo "输入第二个数字: "
#  read anotherNum
#  echo "两个数字分别为 $aNum 和 $anotherNum !"
#  return $(($aNum+$anotherNum))
#}
#funWithReturn
#echo "输入的两个数字之和为 $? !"

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
  echo '-1'
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


a=(1 2 3 4)
tmp_arr=`echo ${a[@]}`
#arrslice "${a[*]}" 2



#test_array=(1 2 3 4)
#test_tmp=`echo ${test_array[@]}`
#num=`arritemidx "$test_tmp" 5`
#echo $num
#unset test_array[${num}]
#echo ${test_array[@]}
arr=(`awk '{print $1}' ${1}|tr '\n' ' '|sed s/' '$//g`)
arr_2=(`awk '{print $1}' ${1}|tr '\n' ' '|sed s/' '$//g`)
for SAMPLE1 in ${arr[@]};do
    #SAMPLE1=`echo $line | awk '{print $1}'`
    #echo ${SAMPLE1}
    #echo ${arr_2[@]}
    for SAMPLE2 in ${arr_2[@]};do
         arr2_tmp=`echo ${arr_2[@]}`
         num=`arritemidx "$arr2_tmp" ${SAMPLE1}`
         #echo $num
         arr_2=(`arrslice "${arr_2[*]}" $num`)
         #if test $num != '-1'; then
         #unset arr_2[${num}]
         #fi
         if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
         echo ${SAMPLE1}_${SAMPLE2}
         fi
    done
done

