#!/usr/bin/env python

"""
@author:yangzz
@file:data_wash.py
@time:2018/10/16
"""

from __future__ import division # very important !!!!
import os
import gzip
import sys
import re

f = open("/data/user/yangzz/mapping/5xresult/combine_3", "r")
lines = f.readlines()
f.close()
dic = {"mid":"1", "high":"2", "low":"0", "NA":"3"}
par_lis = []
for line in lines:
   col = line.strip().split("\t")
   par_lis.append("\t".join([col[0], col[3], dic[col[4]]]))
lis1 = []
lis2 = []
lis3 = []
curr_chr1 = ""
curr_ind1 = -1
curr_chr2 = ""
curr_ind2 = -1
curr_chr3 = ""
curr_ind3 = -1
for line in par_lis:
    col = line.split("\t")
    if re.match('5181',col[1]):
	if col[0] == curr_chr1:
            lis1[curr_ind1].append(col[2])
    	else:
            curr_chr1 = col[0]
            curr_ind1 += 1
            lis1.append([])
            lis1[curr_ind1].append(col[2])
	#
    elif re.match('lx987',col[1]):
        if col[0] == curr_chr2:
            lis2[curr_ind2].append(col[2])
        else:
            curr_chr2 = col[0]
            curr_ind2 += 1
            lis2.append([])
            lis2[curr_ind2].append(col[2])
	#
    else:
        if col[0] == curr_chr3:
            lis3[curr_ind3].append(col[2])
        else:
            curr_chr3 = col[0]
            curr_ind3 += 1
            lis3.append([])
            lis3[curr_ind3].append(col[2])

print(len(lis1))
print(len(lis2))
print(len(lis3))
outpath = "/data/user/yangzz/mapping/5xresult/"
for j in range(0,7):
    for i in range(0,3):
        lett = ["A", "B", "D"]
        filename = outpath + "5181/"  + "chr" + str(j+1) + lett[i]
        f = open(filename, "w")
        for k in lis1[j*3+i]:
            f.write(k+'\n')
        f.close()

        filename = outpath + "lx987/"  + "chr" + str(j+1) + lett[i]
        f = open(filename, "w")
        for k in lis2[j*3+i]:
            f.write(k+'\n')
        f.close()

        filename = outpath + "3097/"  + "chr" + str(j+1) + lett[i]
        f = open(filename, "w")
        for k in lis3[j*3+i]:
            f.write(k+'\n')
        f.close()


#for i in range(0,3):
#    lett = ["A", "B", "D"]
#    filename = "/data/user/yangzz/mapping/test/chrx" + lett[i]
#    f = open(filename, "w")
#    for j in range(0,7):
#        for k in lis[j*3 + i]:
#            f.write(k)
#    f.close()
