#!/bin/bash
a=$1
b=$2
echo $a	$b


#for i in CN_mcl_dir/chr*.*mcl;do
#num=`cat $i | wc -l`
#if [ $num == 0 ];then
#echo $i
#fi
#done

if false;then
for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
(n=$(wc -l < ~/mapping/fieldergenomecompare/statistic/GSR_MCL/${CHR}.1M.bed)
echo $n
#if true;then
for i in `eval echo {1..$n}`;do
echo $i >> ${CHR}_test
done
) &
done
fi

