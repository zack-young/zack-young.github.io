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

dic1 = {}
dic2 = {}
dic3 = {}
dic4 = {}
out1 = open('DP_lx99','w')
out2 = open('GQ_lx99','w')
out3 = open('DP_jm22','w')
out4 = open('GQ_jm22','w')
with open("combine_GQ_lx99", "r") as mydata:
    CNVT = mydata.readlines()
    for line in CNVT:
        item = line.strip().split("\t")
        dic1[item[0]] = dic1.setdefault(item[0], 0)+1
        dic2[item[1]] = dic2.setdefault(item[1], 0)+1
for key,value in dic1.items():
    out1.write('{key}\t{value}\n'.format(key=key, value=value))
out1.close()
for key,value in dic2.items():
    out2.write('{key}\t{value}\n'.format(key=key, value=value))
out2.close()
#
with open("combine_GQ_jm22", "r") as mydata:
    CNVT = mydata.readlines()
    for line in CNVT:
        item = line.strip().split("\t")
        dic3[item[0]] = dic3.setdefault(item[0], 0)+1
        dic4[item[1]] = dic4.setdefault(item[1], 0)+1
for key,value in dic3.items():
    out3.write('{key}\t{value}\n'.format(key=key, value=value))
out3.close()
for key,value in dic4.items():
    out4.write('{key}\t{value}\n'.format(key=key, value=value))
out4.close()


#print(len(dic1["test"]))
#print(dic1["test"])

#            print(ch)
#         print(i,dic1["test"][int(i)][n],dic1["test"][i][n+1])
