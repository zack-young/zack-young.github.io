#!/usr/bin/env python
import re
import sys
#print(1)
CNV1 = open(sys.argv[2],'w')
CNV2 = open(sys.argv[3],'w')
with open(sys.argv[1],'r') as mydata:
    a = 1
    b = 0
    for line in mydata:
        b += 1000000
        Line = line.strip().split("\t")
        priline = line.strip('\n')
        if float(line[0])<0.6 and float(Line[1])<0.6:
            CNV2.write(str(a)+'\t'+str(b)+'\t'+priline+'\t'+'ownCNV'+'\n')
        elif abs(float(Line[0])-float(Line[1]))>0.05:
            CNV1.write(str(a)+'\t'+str(b)+'\t'+priline+'\t'+'unmatchCNV'+'\n')
        #print(line)
        a += 1000000
CNV1.close()
CNV2.close()
