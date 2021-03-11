#!/usr/bin/env python
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

def filter_use(infile, chromosome, filt, compare):
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
    if compare == 'on':
        countl_1 = 0 
        countm_1 = 0
        counth_1 = 0
        countl_2 = 0 
        countm_2 = 0
        counth_2 = 0
        countl_3 = 0 
        countm_3 = 0
        counth_3 = 0
        countl_4 = 0 
        countm_4 = 0
        counth_4 = 0

        for line in IN:
            tokens = line.strip().split("\t") 
            jm22 = str(tokens[3])
            lx99 = str(tokens[4])
            both = str(tokens[5])
            #print(line)
            if jm22 == 'high': 
                
                if lx99 == jm22:
                    if both == 'low_com':
                        countl_1 += 1
                    if both == 'mid_com':
                        countm_1 += 1
                    if both == 'high_com':
                        counth_1 += 1
                else:
                    if both == 'low_com':
                        countl_2 += 1
                    if both == 'mid_com':
                        countm_2 += 1
                    if both == 'high_com':
                        counth_2 += 1                     
            if jm22 == 'low': 
                if lx99 == jm22:
                    if both == 'low_com':
                        countl_3 += 1
                    if both == 'mid_com':
                        countm_3 += 1
                    if both == 'high_com':
                        counth_3 += 1
                        #print(line)
                else:
                    if both == 'low_com':
                        countl_4 += 1
                    if both == 'mid_com':
                        countm_4 += 1
                    if both == 'high_com':
                        counth_4 += 1  
        print("\t".join(['h_h_l', 'h_h_m', 'h_h_h', 
                         'h_l_l', 'h_l_m', 'h_l_h',
                         'l_l_l', 'l_l_m', 'l_l_h',
                         'l_h_l', 'l_h_m', 'l_h_h']))
        print(countl_1, countm_1, counth_1,
              countl_2, countm_2, counth_2,
              countl_3, countm_3, counth_3,
              countl_4, countm_4, counth_4, sep='\t')

    if filt == 'on':
        print("\t".join(['CHR', 'start', 'end', 'alt_ratio_level_jm22', 'alt_ratio_level_lx99']))
        list1 = Makelis(chromosome)
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        start = item[1]
        level = item[3]
    #    IN.seek(0, os.SEEK_SET)
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

    if infile:
        IN.close()

from optparse import OptionParser
# ======================================
def main():
    usage = '''Usage: %prog [-i <input>] [-b 1000] [-o <output>] 
       [-p "jm22_snp_indel_analysis/snp,lx99_snp_indel_analysis/snp"] [-s "_snp_altLEVEL,_snp_altLEVEL"] [-n 2]
            awk '{arr[$4"\t"$5"\t"$6] = arr[$4"\t"$5"\t"$6] + 1}END{for (a in arr) print a"\t"arr[a]}' combine_solo_plus_1.1'''
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
    parser.add_option("--compare", dest="compare",
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
    filter_use(infile=options.infile, chromosome=options.chromosome, filt=options.filt, compare=options.compare)
    #
#
if __name__ == "__main__":
    main()

