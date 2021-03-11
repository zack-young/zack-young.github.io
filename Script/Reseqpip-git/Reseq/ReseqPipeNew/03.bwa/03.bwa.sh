#!/usr/bin/env bash
set -euxo pipefail
# Guo, Weilong; guoweilong@126.com; 2017-10-24

declare -i numb
declare -i numb_MATCH
numb=0
numb_MATCH=0
#

## when files is added later, this number is wrong
for F in ../02.*/$1/??_R1.clean.fq.gz; do
  numb+=1
  ID=`basename $F | cut -c1-2`
  echo $ID
  #wait_yhq
  set +e
  mkdir $1
  set -e
  cd $1
  sh ../bwa.sh ${ID} $1 > ${ID}.o 2 > ${ID}.e
  cd ..
done
#
while [ $numb_MATCH -lt ${numb} ];do
  numb_MATCH=0
  cd $1
  for F in *.o; do
    if cat ${F} | grep -q "Finish"; then
      numb_MATCH+=1
    fi
  done
  cd ..
  echo "Parts total:"${numb}
  echo "Parts done"${numb_MATCH}
  sleep 600
done
#
echo "["`date +%Y-%m-%d,%H:%M:%S`"] All Done"
#rm ../02.*/*.gz
#
#cd ../04.*/
#sh ./04.*.sh
