#!/usr/bin/env python
# -*- coding:utf-8 -*-
from __future__ import division # very important !!!!
import os
import gzip
import sys
import re


def countcolor(infile, s_col=1, chr_col=2):  # 参数顺序要与parser顺序一致！！！
    try:
        if infile :
            if infile.endswith(".gz") :
                IN = gzip.open(infile, 'rb')
            else :
                IN = open(infile, 'r')
            #
        else :
            IN = sys.stdin
        #
    except IOError:
        print("\n[Error]:\n\t File cannot be open: %s" %infile)
        exit(-1)
    #
    # init

    chr_cur = None
    count1 = 0
    count2 = 0
    i = 0
    m = 0
    level = []
    belong = []
    for line in IN:
        tokens = line.strip().split("\t")
        chr = tokens[chr_col - 1]
        if re.match('^chr[1-7][ABD]$', chr):
            level.append(line)
        if re.match('^chr[1-7][ABD].1$', chr):
            belong.append(line)
    for tmp1 in level:
        tokens = tmp1.strip().split("\t")
        gs = tokens[s_col - 1]
        chr = tokens[chr_col - 1]
        if chr == chr_cur:
	    if i == 0:
            	lev = gs
            if i > 0:
                if gs != lev:
                    count1 += 1
		lev = gs
            i += 1
        else:
            if chr_cur is not None:
                print("\t".join([chr_cur, "%s" % count1, "lev"]))
		count1 = 0
		lev = gs
            chr_cur = chr
    print("\t".join([chr_cur, "%s" % count1, "lev"]))
    chr_cur = None
    for tmp2 in belong:
        tokens = tmp2.strip().split("\t")
        gs = tokens[s_col - 1]
        chr = tokens[chr_col - 1]        
        if chr == chr_cur:
	    if m == 0:
            	bel = gs
            if m > 0:
                if gs != bel:
                    count2 += 1
		bel = gs
            m += 1
        else:
            if chr_cur is not None:
                print("\t".join([chr_cur, "%s" % count2, "bel"]))
		count2 = 0
		bel = gs
            chr_cur = chr
    print("\t".join([chr_cur, "%s" % count2, "bel"]))
    if infile:
        IN.close()





from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog [-i <input>] [-b 1000] [-o <output>]\n" \
            "Author : Guo, Weilong; guoweilong@126.com; 2018-06-07\n" \
            "Description: Identify genotype bins orientation from parents.\n" \
            "Input format:\n" \
            "   chr pos geno_p1 geno_p2 geno_son\n" \
            "Output format:\n" \
            "   chr start end total eqp1 eqp2 orient\n" \
            "   chr1A   10001 20000 20  18  2   1\n"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-s", dest="s_col",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 5) # very important!!!
    parser.add_option("-c", dest="chr_col",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    #
    (options, args) = parser.parse_args()
    #
    if (options.outfile is not None) :
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')
        #
    #
    countcolor(options.infile, options.s_col, options.chr_col)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
