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

def GenoPhase(infile, bin_size, highun, lowun, hetedivision, undefinedivision,
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
    count1 = 0
    count2 = 0
    count3 = 0
    count = 0
    #
    l3097 = []
    llx987 = []
    l5181 = []
    for j in range(0,7):
        for i in range(0,3):
            lett = ["A", "B", "D"] 
            l3097.append([])
            l5181.append([])
            llx987.append([])
            f = open("/data/user/yangzz/mapping/5xresult/3097/chr" + str(j+1) + lett[i] + ".1M.norm", "r")
            CNVT = f.readlines()
            for line in CNVT:
                 item = line.strip()
                 l3097[i+j*3].append(item)
            f.close()
            f = open("/data/user/yangzz/mapping/5xresult/lx987/chr" +  str(j+1) + lett[i] + ".1M.norm", "r")
            CNVT = f.readlines()
            for line in CNVT:
                 item = line.strip()
                 llx987[i+j*3].append(item)
            f.close()
            f = open("/data/user/yangzz/mapping/5xresult/5181/chr" +  str(j+1) + lett[i] + ".1M.norm", "r")
            CNVT = f.readlines()
            for line in CNVT:
                 item = line.strip()
                 l5181[i+j*3].append(item)
            f.close()

    for line in IN:
        tokens = line.strip().split("\t")
        chr = tokens[chr_col - 1]
        pos = int(tokens[pos_col - 1])
        g1 = tokens[p1_col - 1]
        g2 = tokens[p2_col - 1]
        gs = tokens[s_col - 1]
        chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
        b = chr_cur[0]
        dic = {"A":"0","B":"1","D":"2"}
        num = int(b[3])-1
        ch = int(dic[b[4]])
        # log((z+1)/(x+y+1)) > 0.1 0.2 0.25
        while pos > end:
            #ratio2 = (undefinedson + 1)/(count + 1)
            #if ratio2 > undefinedivision:
            #    print("\t".join([b, "%s" % start, "%s" % end, "%s" % count, "NA"]))
            #else:
            #    print("\t".join([b, "%s" % start, "%s" % end, "%s" % count, "defined"]))
#	    print("\t".join([b, "%s" % start, "%s" % end, "%s" % hete, "%s" % father, "%s" % mother]))
            count += 1
            if float(l3097[(num)*3+ch][count]) <= 0.6 or float(llx987[(num)*3+ch][count]) <= 0.6 or float(l5181[(num)*3+ch][count]) <= 0.6:
                print("\t".join([b, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2,
                  "%s" % count3,"%s" % count, l3097[(num)*3+ch][count], llx987[(num)*3+ch][count], l5181[(num)*3+ch][count]]))
            start += bin_size
            end += bin_size
            count1 = 0
            count2 = 0
            count3 = 0
            #
            # compare genotype
        if IsGeno(g1) and IsGeno(g2) and IsGeno(gs):                   #continue 会跳出整个循环！！！！
            if g1 == g2:                                               #pass 不起作用 仅是起到站位作用
                pass
            else:
                if IsHomo(g2) and IsHomo(g1):
                    count1 += 2
                else:
                    count1 += 1                       
            #
            if g1 == gs:
                pass
            else:
                if IsHomo(gs) and IsHomo(g1):
                    count2 += 2
                else:
                    count2 += 1
                #
            if g2 == gs:
                pass
            else:
                if IsHomo(g2) and IsHomo(gs):
                    count3 += 2
                else:
                    count3 += 1
        
      #
        #
#    print("\t".join([b, "%s" % start, "%s" % end, "%s" % hete, "%s" % father, "%s" % mother]))
    count += 1
   # print("\t".join([b, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2, "%s" % count3]))
    if float(l3097[(num)*3+ch][count]) <= 0.6 or float(llx987[(num)*3+ch][count]) <= 0.6 or float(l5181[(num)*3+ch][count]) <= 0.6:
        print("\t".join([b, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2,
              "%s" % count3,"%s" % count, l3097[(num)*3+ch][count], llx987[(num)*3+ch][count], l5181[(num)*3+ch][count]]))

    #

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
              int(options.p1_col), int(options.p2_col), int(options.s_col))
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
