#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@author:Levy
@file:finalgeno.py
@time:2018/7/2417:35
"""

from __future__ import division # very important !!!!
from collections import defaultdict
import os
import gzip
import sys
import re
import math
import numpy as np

#awk '$4 > 2 || $4 <1.15 { print $0; }' combine_density > combine_density_threshold

dic1 = defaultdict(list)
f = open("combine_density_threshold", "r")
CNVT = f.readlines()
for j in range(0, 7):
#    print(j)
    for i in range(0, 3):
        dic1["test"].append([])
for line in CNVT:
    item = line.strip().split("\t")
    #print(item[0])
    for j in range(0, 7):
        for i in range(0, 3):
            lett = ["A", "B", "D"]
        #lis = []
        #     lis.append(line)
        # for line in lis:
            if item[0] == "chr"+str(j+1)+lett[i]:
               # print(i+j*3)
                dic1["test"][i+j*3].append(item[1])
                dic1["test"][i+j*3].append(item[2])
        f.close()
#print(len(dic1["test"]))
#print(dic1["test"])

GTF = open("IWGSCv1p1.gtf", "r")
gtf = GTF.readlines()
dic2 = {"A":"0","B":"1","D":"2"}
for line in gtf:
    item = line.strip().split("\t")
    for i in range(len(dic1["test"])):
#    print(len(dic1["test"][int(i)]))     
        ch = int(dic2[item[0][4]])
        if str(i) == 3*item[0][3]+ch:
#            print(ch)
            for n in (np.arange(0,len(dic1["test"][i]),2)):
                if           
#         print(i,dic1["test"][int(i)][n],dic1["test"][i][n+1])
