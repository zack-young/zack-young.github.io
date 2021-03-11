#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from __future__ import division  # very important !!!!
from collections import defaultdict
import os
import gzip
import sys
import re
import math

def Makelis(chromosome):
    global list1
    list1 = []
    f = open(chromosome, "r")
    mask = f.readlines()
    for line in mask:
        item = line.strip("").split("\t")
        list1.append(item[1])
    f.close
    return(list1)

def filter_use(infile, chromosome, filt, ann):
    try:
        if infile:
            if infile.endswith(".gz"):
                IN = gzip.open(infile, 'rb')
            else:
                IN = open(infile, 'r')
            #
        else:
            IN = sys.stdin
        #
    except IOError:
        print("\n[Error]:\n\t File cannot be open: %s" % infile)
        exit(-1)
    if filt == 'on':
        print("\t".join(['CHR', 'start', 'end', 'alt_ratio_level_jm22', 'alt_ratio_level_lx99']))
        for line in IN:
            line_use = line.strip()
            tokens = line_use.split("\t")
            if start in list1:
                pass
            else:
                if start == tokens[1]:
                    print(line_use+'\t'+level)
            start = tokens[1]
            level = tokens[3]

    if ann == 'on':
        print('CHR', 'POS', 'REF', 'ALT','ANN','GENE',
        '954072', 'AD', 'GQ', 'AK58', 'AD', 'GQ', 'fugou583', 'AD', 'GQ',
        'ji200040919', 'AD', 'GQ', 'jimai229', 'AD', 'GQ', 'jimai44', 'AD', 'GQ',
        'jinan17', 'AD', 'GQ', 's742', 'AD', 'GQ', sep='\t')
        for line in IN:
            line_use = line.strip()
            tokens = line.strip().split("\t")
            chr = tokens[0]
            pos = int(tokens[1])
            ref = tokens[2]
            alt = tokens[3]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            # dic2 = {"A":"0","B":"1","D":"2"}
            # num = int(b[3])-1
            # ch = int(dic2[b[4]])
            # count = int(pos/bin_size)
            # len_dic = defaultdict(list)

            ANN = tokens[4].split(",")
            #print(tokens[s1_col - 1],lis1[1],pos,sep='\t')
            #

            lst1=[]
            lst2=[]
            if tokens[4] != '.':
                for ann in ANN:
                    item = ann.split("|")
                    FUN = item[1].strip() #.split("&")
                    lst1.append([FUN, item[2], item[4]])
                [lst2.append(i) for i in lst1 if not i in lst2]
                    #lst2.sort(key=FUN.index)
                for fun in lst2:
                    #print(b,pos,ref,alt,GT1,GT2,"\t".join(item[0:3]),item[4],sep='\t')
                    print(b, pos, ref, alt, fun[0], fun[2], '\t'.join(str(i) for i in tokens[5:]), sep='\t')

    if infile:
        IN.close()

from optparse import OptionParser
# ======================================
def main():
    usage = '''Usage: %prog [-i <input>] [-b 1000] [-o <output>] 
       [-p "jm22_snp_indel_analysis/snp,lx99_snp_indel_analysis/snp"] [-s "_snp_altLEVEL,_snp_altLEVEL"] [-n 2]'''
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                      help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("--chromosome", dest="chromosome",
                      help="column id for offSpring [default: %default]", metavar="STRING",
                      default='1B')
    parser.add_option("--filt", dest="filt",
                      help="column id for offSpring [default: %default]", metavar="STRING",
                      default='off')
    parser.add_option("--ann", dest="ann",
                      help="column id for offSpring [default: %default]", metavar="STRING",
                      default='off')
    parser.add_option("-o", dest="outfile",
                      help="Output file, use STDOUT if omit; "
                      "gzip file end with \".gz\".", metavar="FILE")
    #
    (options, args) = parser.parse_args()
    #
    if (options.outfile is not None) :  #??????
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else:
            sys.stdout = open(options.outfile, 'w')
        #
    #
    filter_use(infile=options.infile, chromosome=options.chromosome, filt=options.filt, ann=options.ann)
    #
#
if __name__ == "__main__":
    main()


