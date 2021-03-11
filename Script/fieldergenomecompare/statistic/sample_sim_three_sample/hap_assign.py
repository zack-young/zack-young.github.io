#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import sys
#~/mapping/fieldergenomecompare/statistic/GSR_MCL/hap_number_count.py

def Entropy ( infile="inflie.txt", file_size=0, sample_lis="sample.txt",CHR="chr1A",pos=1 , pop=10,outfile="outfile.txt"):
    with open(infile) as f:
        GSR_mat=f.readlines()
    num_20=0
    num_10=0
    num_5=0
    num_count_2=0
    if file_size != 0:
        for i in GSR_mat:
            tokens=i.strip().split("\t")
            item_len=len(tokens)
            if item_len > pop*0.2:
                num_20 += 1
            if item_len > pop*0.1 and item_len < pop*0.05 :
                num_10 += 1
            if item_len > pop*0.05 and item_len > 2:
                num_5 += 1
            if item_len <= 2:
                num_count_2 += 1            
    print("\t".join([CHR, pos,"%s" % num_20,"%s" % num_10,"%s" % num_5,"%s" % num_count_2 ]))

from optparse import OptionParser
import math
import numpy as np

# ===========================================
def main():
    parser = OptionParser()
    parser.add_option("-i", dest="infile")
    parser.add_option("-s", dest="sample_lis", default="sample")
    parser.add_option("-f", dest="file_size", default="0")
    parser.add_option("-c", dest="CHR", default="chr1A")
    parser.add_option("-p", dest="pos", default="1")
    parser.add_option("-n", dest="pop", default="1")
    parser.add_option("-o", dest="outfile")
    (options, args) = parser.parse_args()
    #
    Entropy (infile=options.infile, file_size=options.file_size, pop=int(options.pop), sample_lis=options.sample_lis, CHR=options.CHR, pos=options.pos ,outfile=options.outfile)

# ===========================================
if __name__ == "__main__":
    main()
