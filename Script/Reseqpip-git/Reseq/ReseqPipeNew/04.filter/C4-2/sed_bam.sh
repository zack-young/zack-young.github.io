#!/bin/bash
for a in ??.filter.bam;do
#    for num in {1..7}; do
#        for i in {A,B,D};do 
    samtools view -h ${a} |sed "s/_part1/.1/g"|sed "s/_part2/.2/g"|sed "s/Chr/chr/g"|samtools view -Sbh - > ${a}.sed &
        
        #echo "chr${num}${i}part1" >> test.txt
        #echo "chr${num}${i}part2" >> test.txt
#        done
#    done
done
