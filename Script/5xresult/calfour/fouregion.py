#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@author:yangzz
@file:fouregion.py
@time:2018/12/4 17:35
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

def Fouregion(input):
    # open the infile
    try:
        if input:
            IN = open("/data/user/yangzz/mapping/5xresult/"+ input +".raw.bcf.filter", 'r')
            #
        else:
            IN = sys.stdin
        #
    except IOError:
        print("\n[Error]:\n\t File cannot be open: %s" % input)
        exit(-1)
    #
# f = open("/data/user/yangzz/mapping/5xresult/chr1A.fouregion", "r")
# FOUR = f.readlines()
# f.close()
# IN = open("/data/user/yangzz/mapping/5xresult/chr1A.raw.bcf.filter", 'r')
# for line in FOUR:
#     four = line.strip().split("\t")
#     for item in IN:
#         item = item.strip()  # 去除打印时的多余空行！！！
#         vcf = item.strip().split("\t")
#         if vcf[1] >= four[1] and vcf[1] <= four[2]:
#             print(item)
# IN.close()
    f = open("/data/user/yangzz/mapping/5xresult/" + input + ".fouregion", "r")
    FOUR = f.readlines()
    f.close()
    #
    for line in FOUR:
        four = line.strip().split("\t")
        for item in IN:
            item = item.strip()
            vcf = item.strip().split("\t")
            a = int(four[1])
            b = int(four[2])
            c = int(vcf[1])
            if c >= a and c <= b:
                print(item)
            if c > b:
                break  # python以指针方式读取文件，可以用列表方式改善
    #
    if input:
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
    parser.add_option("-i", dest="input",
                  help="Input paramenter, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="STRING")
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
    Fouregion(options.input)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#


# f = open("/data/user/yangzz/mapping/5xresult/chr2A.test_fouregion", "r")
# FOUR = f.readlines()
# f.close()
# IN = open("/data/user/yangzz/mapping/5xresult/chr2A.raw.bcf.filter", 'r')
# for line in FOUR:
#     four = line.strip().split("\t")
#     for item in IN:
#         item = item.strip()  # 去除打印时的多余空行！！！
#         vcf = item.strip().split("\t")
#         a = int(four[1])
#         b = int(four[2])
#         c = int(vcf[1])
#         d = vcf[2]
#         if c >= a and c <= b:
#             print(item)
#         if c > b :
#             break
# #        print("\t".join(["%s" %a,"%s" % b,"%s" % c,"%s" % d]))
# IN.close()