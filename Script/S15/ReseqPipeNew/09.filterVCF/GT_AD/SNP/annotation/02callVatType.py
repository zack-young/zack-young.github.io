#!/usr/bin/env python

import sys
data = open(sys.argv[1])
out = open(sys.argv[2],'w')
out.write('Chr'+'\t'+'pos'+'\t'+'ref'+'\t'+'alt'+'\t'+'ref1_dep'+'\t'+'alt1_dep'+'\t'+'ref2_dep'+'\t'+'alt2_dep'+'\t'+'Var_Type'+'\t'+'sig'+'\t'+'gene'+'\t'+'transcript'+'\t'+'base_Var'+'\t'+'AA_Var'+'\n')
line = data.readlines()
for i in line:
        i = i.split('\t')
	chr1 = i[0].strip()
	pos = i[1].strip()
        info = i[8].strip().split(';')
        info1 = info[-1]
	if info1.startswith('NMD'):
		info2 = info[-3].split(',')
		length = len(info2)
		for j in range(length):
			info3 = info2[j].split('|')
			Type = info3[1].strip()
			sig = info3[2].strip()
			gene = info3[4].strip()
			scr = info3[6].strip()
			AA = info3[9].strip()
			BB = info3[10].strip()
			out.write(chr1+'\t'+str(pos)+'\t'+i[2].strip()+'\t'+i[3].strip()+'\t'+i[4].strip()+'\t'+i[5].strip()+'\t'+i[6].strip()+'\t'+i[7].strip()+'\t'+Type+'\t'+sig+'\t'+gene+'\t'+scr+'\t'+AA+'\t'+BB+'\n')
	elif info1.startswith('LOF'):
		info2 = info[-2].split(',')
		length = len(info2)
                for j in range(length):
                        info3 = info2[j].split('|')
			Type = info3[1].strip()
                        gene = info3[4].strip()
			sig = info3[2].strip()
			scr = info3[6].strip()
			AA =info3[9].strip()
			BB = info3[10].strip()
                        out.write(chr1+'\t'+str(pos)+'\t'+i[2].strip()+'\t'+i[3].strip()+'\t'+i[4].strip()+'\t'+i[5].strip()+'\t'+i[6].strip()+'\t'+i[7].strip()+'\t'+Type+'\t'+sig+'\t'+gene+'\t'+scr+'\t'+AA+'\t'+BB+'\n')
	elif info1.startswith('SOR'):
		continue
	else:
		info4 = info1.split(',')
		length = len(info4)
        	for j in range(length):
                	info5 = info4[j].split('|')
			Type = info5[1].strip()
                	gene = info5[4].strip()
			sig = info5[2].strip()
			scr = info5[6].strip()
			AA = info5[9].strip()
			BB = info5[10].strip()
        		out.write(chr1+'\t'+str(pos)+'\t'+i[2].strip()+'\t'+i[3].strip()+'\t'+i[4].strip()+'\t'+i[5].strip()+'\t'+i[6].strip()+'\t'+i[7].strip()+'\t'+Type+'\t'+sig+'\t'+gene+'\t'+scr+'\t'+AA+'\t'+BB+'\n')

