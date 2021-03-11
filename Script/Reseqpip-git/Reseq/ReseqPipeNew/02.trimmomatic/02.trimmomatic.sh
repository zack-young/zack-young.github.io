#!/usr/bin/env sh
# Guo, Weilong; guoweilong@126.com; 2017-10-24

if [[ $(wc -l ../01.split/ID.list) < 25 ]];then
    yhbatch -n 1 -e trim.e -o trim.o ./trimmomatic.sh ../01.split/ID.list
else
    yhbatch -n 1 -e trim1.e -o trim1.o ./trimmomatic.sh ../01.split/ID.list.part1
    yhbatch -n 1 -e trim2.e -o trim2.o ./trimmomatic.sh ../01.split/ID.list.part2
