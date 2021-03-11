#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import sys
def Entropy ( infile="inflie.txt", file_size=0, sample_lis="sample.txt", outfile="outfile.txt"):
    with open(infile) as f:
        GSR_mat=f.readlines()
    with open(sample_lis) as f:
        sample=f.readlines()
    sample_list = [l.strip().split("\t")[0] for l in sample]
    num=1
    if file_size != 0:
        for i in GSR_mat:
            tokens=i.strip().split("\t")
            for item in tokens:
                print(item+"\t"+"type"+str(num))
                sample_list.remove(item)
            num+=1
    for i in sample_list:
        print(i+"\t"+"other")
        num+=1
from optparse import OptionParser
import math
import numpy as np

# ===========================================
def main():
    parser = OptionParser()
    parser.add_option("-i", dest="infile")
    parser.add_option("-s", dest="sample_lis", default="sample")
    parser.add_option("-f", dest="file_size", default="0")
    parser.add_option("-o", dest="outfile")
    (options, args) = parser.parse_args()
    #
    Entropy (infile=options.infile, file_size=options.file_size,sample_lis=options.sample_lis,  outfile=options.outfile)

# ===========================================
if __name__ == "__main__":
    main()
