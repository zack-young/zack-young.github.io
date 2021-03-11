#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@author:Levy
@file:finalgeno.py
@time:2018/7/2417:35
"""

from __future__ import division # very important !!!!
import os
import gzip
import sys
import re
import math

def error(msg):
    print(sys.stderr, 'ERROR: %s' % msg)
    exit(1)
#
def IsGeno(geno):
    if len(geno)== 3 and geno[1] =='/' and geno!= "./.":
        return 1
    else:
        return 0
    #
#
def IsHete(geno):
    if geno[0] == geno[2]:
        return 0
    else:
        return 1
    #
#
# TODO: 增加更改参数功能
def IsHomo(geno):  #判断纯合
    if geno[0] == geno[2]:
        return 1
    else:
        return 0

def  Isconsistent(g1, g2):   #判断父母本差异密度
    if g1 == g2: #and g2 == gs :
        return 1
    #
    # if IsGeno(g1) and IsGeno(g2) and IsGeno(gs) :
    #     return 1
    #   
    else :
        return 0
    #
#

def OutputLinehomohigh(chr_cur, start, end, count, count1, count2) :
    if count1 == count2 :
        state = "NA"
    elif count1 > count2 :
        state = "high987"
    else :
        state = "high3097"
    #
    print("\t".join([chr_cur, "%s" % start, "%s" % end,
                     "%s" % count, "%s" % count1, "%s" % count2, state]) )
    #
#

def OutputLinehomomid(chr_cur, start, end, count, count1, count2) :
    if count1 == count2:
        state = "NA"
    elif count1 > count2:
        state = "mid987"
    else:
        state = "mid3097"
    #
    print("\t".join([chr_cur, "%s" % start, "%s" % end,
                     "%s" % count, "%s" % count1, "%s" % count2, state]) )
    #
#

def GenoPhase (infile, bin_size, highun, lowun, division,
               chr_col=1, pos_col=2, p1_col=3, p2_col=4, s_col=5):
    # open the infile
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
    #
    # init
    start = 1
    end = start + bin_size - 1
    count = 0
    count1 = 0
    count2 = 0
    consistent = 0
    unconsis = 0
    hete = 0
    #
    for line in IN:
        tokens = line.strip().split("\t")
        chr = tokens[chr_col - 1]
        pos = int(tokens[pos_col - 1])
        g1 = tokens[p1_col - 1]
        g2 = tokens[p2_col - 1]
        gs = tokens[s_col - 1]
        chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
        b = chr_cur[0]
        # log((z+1)/(x+y+1)) > 0.1 0.2 0.25
        if IsGeno(g1) and IsGeno(g2) and IsGeno(gs) and IsHomo(g1) and IsHomo(g2):
            while pos > end:
                ratio = math.log((hete + 1)/(count1 + count2 + 1))
                if unconsis >= highun:
                    if ratio > division:
                        print("\t".join([b , "%s" % start, "%s" % end, "%s" % count, "NA", "NA", "highhete"]))
                    else:
                        OutputLinehomohigh(b, start, end, count, count1, count2)
                elif lowun <= unconsis < highun:
                    if ratio > division:
                        print("\t".join([b, "%s" % start, "%s" % end, "%s" % count, "NA", "NA", "midhete"]))
                    else:
                        OutputLinehomomid(b, start, end, count, count1, count2)
                else:
                    print("\t".join([b, "%s" % start, "%s" % end, "%s" % count, "NA", "NA", "low"]))
                start += bin_size
                end += bin_size
                count = 0
                count1 = 0
                count2 = 0
                consistent = 0
                unconsis = 0
                hete = 0
                #
                #
            #
            # compare genotype
            count += 1
            
            if IsHomo(gs):
                if g1 == gs:
                    count1 += 1
                else:
                    count2 += 1
            else:
                hete += 1
            if Isconsistent(g1, g2):
                consistent += 1
            else:
                unconsis += 1
            #
        #
        #

    #
    
    if IsGeno(g1) and IsGeno(g2) and IsGeno(gs) and IsHomo(g1) and IsHomo(g2):
        if IsHomo(gs):
            if g1 == gs:
                count1 += 1
            else:
                count2 += 1
        else:
            hete += 1        
        if Isconsistent(g1, g2):
            consistent += 1
        else:
            unconsis += 1
        ratio = math.log10((hete + 1)/(count1 + count2 + 1))
        if unconsis >= highun:
            if ratio > division:
                print("\t".join([b, "%s" % start, "%s" % end, "%s" % count, "NA", "NA", "highhete"]))
            else:
                OutputLinehomohigh(b, start, end, count, count1, count2)
        elif lowun <= unconsis < highun:
            if ratio > division:
                print("\t".join([b, "%s" % start, "%s" % end, "%s" % count, "NA", "NA", "midhete"]))
            else:
                OutputLinehomomid(b, start, end, count, count1, count2)
        else:
            print("\t".join([b, "%s" % start, "%s" % end, "%s" % count, "NA", "NA", "low"]))


    #
    if infile:
        IN.close()
    #
#



from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog [-i <input>] [-b 1000] [-o <output>]\n" \
            "Author : Yang, zhengzhao; yangzhengzhao@cau.edu.cn; 2018-06-07\n" \
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
    parser.add_option("-b", dest="bin_size",
                  help="size of bin (bp) [default: %default]", metavar="INT",
                      default = 10000)
    parser.add_option("-m", dest="highun",
                  help="minimum of high inconsistent [default: %default]", metavar="INT",
                      default = 2000)
    parser.add_option("-l", dest="lowun",
                  help="maximum of low inconsistent [default: %default]", metavar="INT",
                      default = 100)
    parser.add_option("-d", dest="division",
                  help="division [default: %default]", metavar="INT",
                      default = 0.1)                  
    parser.add_option("-c", dest="chr_col",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-p", dest="pos_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 2)
    parser.add_option("-1", dest="p1_col",
                  help="column id for 1st parent [default: %default]", metavar="INT",
                      default = 5)
    parser.add_option("-2", dest="p2_col",
                  help="column id for 2nd parent [default: %default]", metavar="INT",
                      default = 7)
    parser.add_option("-s", dest="s_col",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 6)
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
    GenoPhase (options.infile, int(options.bin_size),int(options.highun),int(options.lowun),float(options.division),
               options.chr_col, options.pos_col,
               options.p1_col, options.p2_col, options.s_col)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
