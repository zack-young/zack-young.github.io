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
def data_wash(infile, outpath):
    IN = open(infile, "r")
    lines = IN.readlines()
    dic1 = {"lx987":"1", "987":"1", "3097":"2", "hete":"0", "NA":"3"}  #  !!!!!!!!!
    dic2 = {"mid":"1", "high":"2", "low":"0", "NA":"3"}
    par_lis1 = []
    par_lis2 = []
    for line in lines:
       col = line.strip().split("\t")
       if col[0].find(".1") != -1 :
           par_lis1.append("\t".join([col[0], dic1[col[4]]]))
       if len(col[0]) == 5 :
           par_lis2.append("\t".join([col[0], dic2[col[4]]])) 
    lis1 = []
    curr_chr = ""
    curr_ind = -1
#######
    for line in par_lis1:
        col = line.split("\t")
        if col[0] == curr_chr:
            lis1[curr_ind].append(col[1])
        else:
            curr_chr = col[0]
            curr_ind += 1
            lis1.append([])
            lis1[curr_ind].append(col[1])
    len(lis1)
    print(len(lis1))
    for j in range(0,7):   #  (0,3) 0,1,2
        for i in range(0,3):
            lett = ["A", "B", "D"]
            filename = outpath + "belong/" + "chr" + str(j+1) + lett[i]
            f = open(filename, "w")
            for k in lis1[j*3+i]:
                f.write(k+'\n')
            f.close()
#######
    lis2 = []
    curr_chr = ""
    curr_ind = -1
    for line in par_lis2:
        col = line.split("\t")
        if col[0] == curr_chr:
            lis2[curr_ind].append(col[1])
        else:
            curr_chr = col[0]
            curr_ind += 1
            lis2.append([])
            lis2[curr_ind].append(col[1])
    len(lis2)
    print(len(lis2))
    for j in range(0,7):
        for i in range(0,3):
            lett = ["A", "B", "D"]
            filename = outpath + "density/" + "chr" + str(j+1) + lett[i]
            f = open(filename, "w")
            for k in lis2[j*3+i]:
                f.write(k+'\n')
            f.close()
  
    if infile:
	IN.close()
#    for i in range(0,3):
#        lett = ["A", "B", "D"]
#        filename = "/data/user/shinyug/python/HMMinuse/data/5M_parents/1x3/chrx" + lett[i]
#        f = open(filename, "w")
#        for j in range(0,7):
#            for k in lis[j*3 + i]:
#                f.write(k)
#        f.close()
from optparse import OptionParser

def main():
    usage = "Usage: %prog [-i <input>] [-b 1000] [-o <output>]\n"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-o", dest="outpath",
                  help="Output path, use STDOUT if omit; ", metavar="PATH")
    (options, args) = parser.parse_args()

    data_wash(options.infile, options.outpath)

if __name__ == "__main__":
    main()

