#!/usr/bin/env python
from __future__ import print_function
import sys
import re
data = open(sys.argv[1],'r')
#out = open(sys.argv[2],'w')
print('CHR','POS','REF','ALT','lx99','lx99_readdeepth','jm22','jm22_readdeepth','ANN',sep=',')
for item in data:
    if not item.startswith('#'):
        i = item.strip().split('\t')
        chr1 = i[0]
        pos = i[1]
        alt = i[4].split(',')
        if i[7].split(';')[-1].startswith('SOR'):        
            continue
        else:
            ann = i[7].split('ANN=')[1]
        spec_ann = ann.split(',')
        sample1 = i[9].split(':')[0][0]
        sample2 = i[10].split(':')[0][0]
        sample1_1 = i[9].split(':')[0]
        sample2_2 = i[10].split(':')[0]
        DP1 = i[9].split(':')[2]
        DP2 = i[10].split(':')[2]
        num1 = 0
        num2 = 0
        if int(sample1) != 0:
            index1 = int(sample1)-1
            num1 = 1
        if int(sample2) != 0:
            index2 = int(sample2)-1
            num2 = 1
        if num1 == 1:
            alt1 = alt[index1]
        else:
            alt1 = 'nun'
        if num2 == 1:
            alt2 = alt[index2]
        else:
            alt2 = 'num'    
        list1 = []
        for info in spec_ann:
            infobase = info.split('|')
            base = infobase[0]
            if base == alt1 or base == alt2:
                list1.append(info)
        for item in list1:
#        print(chr1, pos, i[4], sample1_1, sample2_2, *list1, sep='\t')
            print(chr1, pos, i[3], i[4], sample1_1,DP1, sample2_2,DP2, item, sep=',')
data.close
