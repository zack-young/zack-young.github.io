#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@author:Levy
@file:finalgeno.py
@time:2018/7/2417:35
"""

from __future__ import division # very important !!!!
from collections import defaultdict
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
def  Isconsistent(g1, g2):   #判断父母本差异密度
    if g1 == g2: #and g2 == gs :
        return 1
    # if IsGeno(g1) and IsGeno(g2) and IsGeno(gs) :
    #     return 1  
    else :
        return 0
    #
#

def OutputLinehomohigh(chr_cur, start, end, count, count1, count2) :
    if count1 == count2 :
        state = "NA"
    elif count1 > count2 :
        state = "lx987"
    else :
        state = "3097"
    #
    print("\t".join([chr_cur, "%s" % start, "%s" % end,
                     "%s" % count, "%s" % count1, "%s" % count2, state]) )
    #
#

def OutputLinehomomid(chr_cur, start, end, count, count1, count2) :
    if count1 == count2:
        state = "NA"
    elif count1 > count2:
        state = "lx987"
    else:
        state = "3097"
    #
    print("\t".join([chr_cur, "%s" % start, "%s" % end,
                     "%s" % count, "%s" % count1, "%s" % count2, state]) )
    #
#
def DefineAD(gt, ad):
    global lis1
    lis1 = []
    item_ad = ad.strip().split(",")
    if IsHete(gt):
        if int(item_ad[0]) == 1:
            gt = "1/1"
        if int(item_ad[1]) == 1:
            gt = "0/0"
    if ad != '.':
        lis1.append(int(item_ad[0])+int(item_ad[1]))
    else:
        lis1.append(0)
    lis1.append(gt)
    return(lis1)
#
def ChangeGQ(gt, gq):
    global lis2
    lis2 = []
    if IsGeno(gt):
        lis2.append(int(gq))
    else:
        lis2.append(0)
    return(lis2)
#
def ChangeGQ_dot(gq):
    global lis3
    lis3 = []
    if gq != '.':
        lis3.append(int(gq))
    else:
        lis3.append(0)
    return(lis3)
#
def Makelis(chromsome):
    global list1
    list1 = []
    #f = open("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/MASK/chr" + chromsome + "_combine_total", "r")
    f = open(chromsome, "r")
    mask = f.readlines()
    for line in mask:
        item = line.strip("").split("\t")
        list1.append(item[1])
    f.close
    return(list1)
#
def Makedic(sample):
    dic1 = defaultdict(list)
    for j in range(0, 7):
        for i in range(0, 3):
            lett = ["A", "B", "D"]
            dic1[sample].append([])
            f = open(sample + "/chr" + str(j+1) + lett[i] + ".1M.norm", "r")
            CNVT = f.readlines()
            f.close()
            #lis = []
            for line in CNVT:
            #     lis.append(line)
            # for line in lis:
                item = line.strip()
                dic1[sample][i+j*3].append(item)
    return(dic1)


def GenoPhase(infile, bin_size, mult, single_DP_GQ, chr_col=1, pos_col=2, 
              ref_col=3, alt_col=4, gt_col=5, p2_col=4, s_col=6,  chromsome='1A'):
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

    # firstline = IN.readline().strip()
    # item = firstline.strip().split("\t")
    # firstpos = int(item[pos_col - 1])
    # if firstpos <= bin_size:
    #     start = 1
    #     end = bin_size
    # else:
    #     s_num = int(firstpos/bin_size)-1
    #     start = (s_num+1) * bin_size
    #     end = start + bin_size
    # IN.seek(0, os.SEEK_SET) 
    if mult == 'on':   # the function of this module is to transform AD to DP
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            #
            len_dic = defaultdict(list)
            list1 = []
            a = 0
            for i in range(int((len(tokens)-4)/3)):
                #print(tokens,'\t', '%s' % i, tokens[4+i*3])
                GT1 = tokens[4+i*3]
                AD1 = tokens[5+i*3]
                item_ad1 = AD1.strip().split(",")
                GQ1 = tokens[6+i*3]
                ChangeGQ(GT1, GQ1)
                GQ1 = lis2[0]
                DefineAD(GT1, AD1)
                DP1 = lis1[0]
                #len_dic[i].append([GT1, DP1, GQ1])
                list1.append(GT1)
                list1.append(DP1)
                list1.append(GQ1)
                if IsGeno(GT1):
                    a += 1
            #
            if a > 0:
                #print("\t".join([b, '%s' % pos, ref, alt, *list1]))
                print(b, pos, ref, alt, *list1, sep='\t')
    if single_DP_GQ == 'on':
        '''
        this module is created to product DP counts and GQ counts in miss homo hete variants
        for determine threshold
        '''
        f = open(chromsome, "r")
        mask = f.readlines()
        f.close
        a = 0
        dic1 = {}
        dic1_1 = {}                
        dic2 = {}
        dic2_2 = {}
        dic3 = {}
        dic3_3 = {}
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            #
            len_dic = defaultdict(list)
            a = 0
            GT = tokens[s_col - 1]
            AD = tokens[s_col]
            GQ = tokens[s_col + 1]
            lis3 = ChangeGQ_dot(GQ)
            GQ = lis3[0]
            lis1 = DefineAD(GT, AD)
            DP = lis1[0]
            GT = lis1[1]
            judge_num = 0
            for ther in mask:
                item = ther.strip("").split("\t")
                if item[1] <= str(pos) <= item[2]:
                    judge_num += 1       
            if judge_num == 0:
                if IsGeno(GT):
                    if IsHete(GT):
                        dic1[DP] = dic1.setdefault(DP, 0)+1
                        dic1_1[GQ] = dic1_1.setdefault(GQ, 0)+1
                    else:
                        dic2[DP] = dic2.setdefault(DP, 0)+1
                        dic2_2[GQ] = dic2_2.setdefault(GQ, 0)+1
                else:
                    dic3[DP] = dic3.setdefault(DP, 0)+1  
                    dic3_3[GQ] = dic3_3.setdefault(GQ, 0)+1      
        for key,value in dic1.items():
            print('{key}\t{value}'.format(key=key, value=value),'hete_DP',sep='\t')
        for key,value in dic1_1.items():
            print('{key}\t{value}'.format(key=key, value=value),'hete_GQ',sep='\t')
        for key,value in dic2.items():
            print('{key}\t{value}'.format(key=key, value=value),'homo_DP',sep='\t')
        for key,value in dic2_2.items():
            print('{key}\t{value}'.format(key=key, value=value),'homo_GQ',sep='\t')
        for key,value in dic3.items():
            print('{key}\t{value}'.format(key=key, value=value),'miss_DP',sep='\t')
        for key,value in dic3_3.items():
            print('{key}\t{value}'.format(key=key, value=value),'miss_GQ',sep='\t')            







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
    parser.add_option("--mult", dest="mult",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--single_DP_GQ", dest="single_DP_GQ",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("-c", dest="chr_col",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-p", dest="pos_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 2)
    parser.add_option("-r", dest="ref_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 3)
    parser.add_option("-a", dest="alt_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 4)
    parser.add_option("-1", dest="gt_col",
                  help="column id for 1st parent [default: %default]", metavar="INT",
                      default = 5)
    parser.add_option("-2", dest="p2_col",
                  help="column id for 2nd parent [default: %default]", metavar="INT",
                      default = 7)
    parser.add_option("-s", dest="s_col",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 6)
    parser.add_option("--chromsome", dest="chromsome",
                  help="column id for offSpring [default: %default]", metavar="STRING",
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
    GenoPhase(options.infile, int(options.bin_size), options.mult, options.single_DP_GQ, int(options.chr_col), int(options.pos_col),
              int(options.ref_col), int(options.alt_col),
              int(options.gt_col), int(options.p2_col), int(options.s_col), options.chromsome)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
