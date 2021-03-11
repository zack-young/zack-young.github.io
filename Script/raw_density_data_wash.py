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

f = open("/data/user/yangzz/mapping/5xresult/5181/threelevel", "r")
lines = f.readlines()
f.close()
dic = {"mid":"1", "high":"2", "low":"0", "NA":"3"}
par_lis = []
for line in lines:
   col = line.strip().split("\t")
   if len(col[0]) == 5 :
       par_lis.append("\t".join([col[0], dic[col[3]]]))
lis = []
curr_chr = ""
curr_ind = -1
for line in par_lis:
    col = line.split("\t")
    if col[0] == curr_chr:
        lis[curr_ind].append(col[1])
    else:
        curr_chr = col[0]
        curr_ind += 1
        lis.append([])
        lis[curr_ind].append(col[1])
len(lis)
for j in range(0,7):
    for i in range(0,3):
        lett = ["A", "B", "D"]
        filename = "/data/user/yangzz/mapping/5xresult/5181/chr" + str(j+1) + lett[i]
        f = open(filename, "w")
        for k in lis[j*3+i]:
            f.write(k+'\n')
        f.close()

for i in range(0,3):
    lett = ["A", "B", "D"]
    filename = "/data/user/yangzz/mapping/test/chrx" + lett[i]
    f = open(filename, "w")
    for j in range(0,7):
        for k in lis[j*3 + i]:
            f.write(k)
    f.close()
