#!/usr/bin/env python
# -*- coding:utf-8 -*-

#dic1 = defaultdict(list)
#for j in range(0, 7):
#    for i in range(0, 3):
#        lett = ["A", "B", "D"]
#        dic1[sample].append([])
#        f = open( sample + "/chr" + str(j+1) + lett[i] + ".1M.norm", "r")
#        CNVT = f.readlines()
#        #lis = []
#        for line in CNVT:
#        #     lis.append(line)
#        # for line in lis:
#            item = line.strip()
#            dic1[sample][i+j*3].append(item)
#        f.close()
f = open("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/MASK/chr" + chromosome + ".mask_filtered_combine", "r")
list1 = []
mask = f.readlines()
#lis = []
for line in CNVT:
#     lis.append(line)
# for line in lis:
    item = line.strip()
    dic1[sample][i+j*3].append(item)
f.close()
