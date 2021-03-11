#!/usr/bin/env python
import sys
data = open(sys.argv[1])
#out = open(sys.argv[2],'w')
da = data.readlines()
dict1 = {}
for i in da:
    if not i.startswith('#'):
        i = i.split('\t')
        ann = i[7].split('ANN')
        print(ann[0]+'\t'+ann[1])
#	if i[5].strip() != 'MODIFIER':
#		gene = i[2]
#		dict1.setdefault(gene,[]).append(i[4].strip()+' '+i[5].strip())
#for key in dict1.keys():
#	info = dict1[key]
#	Len = len(info)
#	info1 = ''
#	for j in range(Len):
#		info1 = info1 + info[j] + ' '
#	out.write(key+'\t'+info1+'\n')
