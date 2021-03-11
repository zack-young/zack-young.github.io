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
def  IsHomo(geno):  #判断纯合
    if geno[0] == geno[2]:
        return 1
    else:
        return 0
    #
#
def  Isunconsistent(g1, g2):   #判断父母本差异密度
    if g1 == g2: #and g2 == gs :
        return 0
    # if IsGeno(g1) and IsGeno(g2) and IsGeno(gs) :
    #     return 1  
    else:
        return 1
    #
#
def Makelis(chromsome):
    global list1
    list1 = []
    f = open("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/MASK/lx99_jm22_mask/chr" + chromsome + ".mask_filtered_combine", "r")
    mask = f.readlines()
    for line in mask:
        item = line.strip("").split("\t")
        list1.append(item[2])
    f.close
    #
#

def compare(b, start, end, hete, father, mother, axis, ayis) :
    if axis > 2.84 and ayis <= 2.89:
        print("\t".join([b, "%s" % start, "%s" % end, "5","%s" % hete, "%s" % father, "%s" % mother]))
    #
    if axis <= 2.84 and ayis > 2.89:
        print("\t".join([b, "%s" % start, "%s" % end, "1","%s" % hete, "%s" % father, "%s" % mother]))
    if axis < 2.84 and ayis < 2.89:
        if ayis > 2.3/3.5*axis + 1:
            print("\t".join([b, "%s" % start, "%s" % end, "2","%s" % hete, "%s" % father, "%s" % mother]))
        elif ayis < 5/3.5*axis - 1.2 :
            print("\t".join([b, "%s" % start, "%s" % end, "4","%s" % hete, "%s" % father, "%s" % mother]))
        else:
            print("\t".join([b, "%s" % start, "%s" % end, "3","%s" % hete, "%s" % father, "%s" % mother]))

def GenoPhase(infile, bin_size, highun, lowun, hetedivision, undefinedivision,
              chr_col=1, pos_col=2, p1_col=3, p2_col=4, s_col=5, chromsome='1A'):
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
    firstline = IN.readline().strip()
    item = firstline.strip().split("\t")
    firstpos = int(item[pos_col - 1])
    if int(firstpos) <= int(bin_size):
        start = 1
        end = bin_size
    else:
        num = int(firstpos/bin_size)-1
        start = (num+1) * bin_size
        end = start + bin_size
    IN.seek(0, os.SEEK_SET)
    count1 = 0
    count2 = 0
    count3 = 0
    hete = 0
    son = 0
    father = 0
    mother = 0
    undefinedpar = 0
    Makelis(chromsome)
#    print(list1)
#    continue
    #
    for line in IN:
        tokens = line.strip().split("\t")
        chr = tokens[chr_col - 1]
        pos = int(tokens[pos_col - 1])
        g1 = tokens[p1_col - 1]
        g2 = tokens[p2_col - 1]
        chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
        b = chr_cur[0]
        while pos > end:
            START = str(start)
            deli = abs(count2-count1)
            ratio = count3/deli
            if START in list1:
                pass
            else:
                if count3 >= 4063:
                    print("\t".join([b, "%s" % start, "%s" % end, "%.4f" % ratio, "%s" % count3, "%s" % count1, "level1"]))
                elif count3 < 974:
                    print("\t".join([b, "%s" % start, "%s" % end, "%.4f" % ratio, "%s" % count3, "%s" % count1, "level4"]))
                elif 974 <= count3 < 2362:
                    print("\t".join([b, "%s" % start, "%s" % end, "%.4f" % ratio, "%s" % count3, "%s" % count1, "level3"]))
                else:
                    print("\t".join([b, "%s" % start, "%s" % end, "%.4f" % ratio, "%s" % count3, "%s" % count1, "level2"]))
            # #
            start += bin_size
            end += bin_size
            count1 = 0
            count2 = 0
            count3 = 0
            hete = 0
            son = 0
            father = 0
            mother = 0
            #
        if IsGeno(g1):
            count1 += 1
        if IsGeno(g2):
            count2 += 1
        if IsGeno(g1) or IsGeno(g1):
            if Isunconsistent(g1, g2):
                count3 += 1
       #
#    print("\t".join([b, "%s" % start, "%s" % end, "%s" % hete, "%s" % father, "%s" % mother]))
#    if father > mother:
#        ratio = father/(mother+1)
#    else:
#        ratio = mother/(father+1)
    deli = abs(count2-count1)
    ratio = count3/deli
    START = str(start)
    if START in list1:
        pass
    else:
        if count3 >= 4063:
            print("\t".join([b, "%s" % start, "%s" % end, "%.4f" % ratio, "%s" % count3, "%s" % count1, "level1"]))
        elif count3 < 974:
            print("\t".join([b, "%s" % start, "%s" % end, "%.4f" % ratio, "%s" % count3, "%s" % count1, "level4"]))
        elif 974 <= count3 < 2362:
            print("\t".join([b, "%s" % start, "%s" % end, "%.4f" % ratio, "%s" % count3, "%s" % count1, "level3"]))
        else:
            print("\t".join([b, "%s" % start, "%s" % end, "%.4f" % ratio, "%s" % count3, "%s" % count1, "level2"]))
#        print("\t".join([b, "%s" % start, "%s" % end, "%s" % hete, "%s" % father, "%s" % mother]))
#    print("\t".join([b, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2, "%s" % count3]))
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
    parser.add_option("-d", dest="hetedivision",
                  help="division [default: %default]", metavar="INT",
                      default = 0.1)
    parser.add_option("-u", dest="undefinedivision",
                  help="division [default: %default]", metavar="INT",
                      default = 0.42)                  
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
    parser.add_option("-n", dest="chromsome",
                  help="chromsome number [default: %default]", metavar="STRING",
                      default = '1A')
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
    GenoPhase(options.infile, int(options.bin_size), int(options.highun), int(options.lowun), 
              float(options.hetedivision), float(options.undefinedivision),
              options.chr_col, options.pos_col,
              int(options.p1_col), int(options.p2_col), int(options.s_col), str(options.chromsome))
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
