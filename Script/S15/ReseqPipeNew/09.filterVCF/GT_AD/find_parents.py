#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@author:Levy
@file:finalgeno.py
@time:2018/7/2417:35
"""

from __future__ import division  # very important !!!!
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
def IsALT(geno):
    if geno[0] != '0' or geno[2] != "0":
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
    else :
        return 1
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
        if '1' in item_ad:
            index = [i for i,x in enumerate(item_ad) if x != '1' and x != '0'] # find all not '1' and '0' item index
            if len(index) == 1:
                gt = str(index[0])+'/'+str(index[0])
            elif len(index) == 0:
                gt = '0/0'
            else:  
                gt = str(index[0])+'/'+str(index[1])
    if ad != '.':
        if gt != './.':
            lis1.append(int(item_ad[int(gt[0])])+int(item_ad[int(gt[2])]))
        else:
            lis1.append(sum(map(int, item_ad)))
    else:
        lis1.append(0)
    lis1.append(gt)
    return(lis1)

def ChangeGQ(gt, gq):
    global lis2
    lis2 = []
    if IsGeno(gt):
        lis2.append(int(gq))
    else:
        lis2.append(0)
    return(lis2)

def ChangeGQ_dot(gq):
    global lis3
    lis3 = []
    if gq != '.':
        lis3.append(int(gq))
    else:
        lis3.append(0)
    return(lis3)

def Makelis(chromsome):
    global list1
    list1 = []
    #f = open("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/MASK/chr" + chromsome + "_combine_total", "r")
    f = open(chromsome, "r")
    mask = f.readlines()
    for line in mask:
        item = line.strip("\n").split("\t")
        list1.append([item[1],item[2], item[3]])
    f.close
    return(list1)

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
def GenoPhase(infile, bin_size, chromosome, pre_define,LEVEL,DP_high, DP_low, GQ_sample,
              chr_col=1, pos_col=2, ref_col=3, alt_col=4, s1_col=5, s2_col=4, s3_col=5):
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
    count1 = 0
    count2 = 0
    count1_1 = 0
    count2_1 = 0
    count3 = 0
    count4 = 0
    count5 = 0
    count6 = 0
    count7 = 0
    count8 = 0
    count9 = 0
    count5_1 = 0
    count6_1 = 0
    count7_1 = 0
    count8_1 = 0
    DP_count = 0
    length = {'chr1A': '594102056', 'chr1B': '689851870', 'chr1D': '495453186',
              'chr2A': '780798557', 'chr2B': '801256715', 'chr2D': '651852609',
              'chr3A': '750843639', 'chr3B': '830829764', 'chr3D': '615552423', 
              'chr4A': '744588157', 'chr4B': '673617499', 'chr4D': '509857067',
              'chr5A': '709773743', 'chr5B': '713149757', 'chr5D': '566080677',
              'chr6A': '618079260', 'chr6B': '720988478', 'chr6D': '473592718',
              'chr7A': '736706236', 'chr7B': '750620385', 'chr7D': '638686055',
              'chrUn': '480980714'}
#    Makedic(sample)
    #
    if pre_define == 'on':   
        f = open(chromosome, "r")
        mask = f.readlines()[1:]
        f.close
        define_list=[]
        for line in mask:
            item = line.strip().split("\t")
            if 'CNV' in str(item[3]) or 'CNV' in str(item[4]):
                if str(item[3]) == 'both_CNV':
                    if str(item[4]) == 'both_CNV':
                        define_list.append(item[0:3]+["undefined_CNV"])
                    else:
                        define_list.append(item[0:3]+["954072_CNV"])
                elif str(item[4]) == 'both_CNV':
                    define_list.append(item[0:3]+["JiNan17_CNV"])
                else:
                    define_list.append(item[0:3]+["undefined_CNV"])             
            else:
                if str(item[3]) == 'high_com':
                    if str(item[4]) == 'low_com' or str(item[4]) == 'mid_com':
                        define_list.append(item[0:3] + ["JiNan17_snp"])
                    else:
                        define_list.append(item[0:3]+["undefined_snp"])
                elif str(item[3]) == 'mid_com':
                    if str(item[4]) == 'low_com':
                        define_list.append(item[0:3] + ["JiNan17_snp"])
                    elif str(item[4]) == 'high_com':
                        define_list.append(item[0:3]+["954072_snp"])
                    else:
                        define_list.append(item[0:3]+["wait_define_mid"])
                else:
                    if str(item[4]) == 'low_com':
                        define_list.append(item[0:3] + ["wait_define_low"])
                    else:
                        define_list.append(item[0:3]+["954072_snp"])
        print("\t".join(['CHR', 'start', 'end', 'diff_homosnp_ratio_level']))
        #for i in range(0,len(define_list)):
        #    print("\t".join(str(i) for i in define_list[i][0:4]))
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        count_same_gt1_gt2 = 0
        count_same_gt2_gt3 = 0
        gt1_gts_same_num = 0
        gt2_gts_same_num = 0 
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        # PART = infile.split("_")
        # part = PART[0]
        # art = part.lstrip('s')
        #level = LEVEL.strip().split(",")
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            #
            GT1 = tokens[s1_col - 1]
            AD1 = tokens[s1_col]
            GQ1 = tokens[s1_col + 1]
            lis3 = ChangeGQ_dot(GQ1)
            GQ1 = lis3[0]
            lis1 = DefineAD(GT1, AD1)
            DP1 = lis1[0]
            GT1 = lis1[1]
            #
            GT2 = tokens[s2_col - 1]
            AD2 = tokens[s2_col]
            GQ2 = tokens[s2_col + 1]
            lis3_1 = ChangeGQ_dot(GQ2)
            GQ2 = lis3_1[0]
            lis1_1 = DefineAD(GT2, AD2)
            DP2 = lis1_1[0]
            GT2 = lis1_1[1]
            #
            GTs = tokens[s3_col - 1]
            ADs = tokens[s3_col]
            GQs = tokens[s3_col + 1]
            lis3_2 = ChangeGQ_dot(GQs)
            GQs = lis3_2[0]
            lis1_2 = DefineAD(GTs, ADs)
            DPs = lis1_2[0]
            GTs = lis1_2[1]
            while pos > end:
                START = str(start)
                judge = 0
                #deli = abs(count2-count1)
                for i in range(0, len(define_list)):
                    if START == define_list[i][1]:
                        if define_list[i][3] != 'wait_define_low':
                            #print("\t".join(str(i) for i in define_list[i][0:4]))
                        #print("break1")
                            judge = 1
                        else:
                            judge = 2
                        #print(len(define_list),START ,define_list[i][1],define_list[i][3],judge)
                if judge == 2:
                    count_same_gt1_gt2 += gt1_gts_same_num
                    count_same_gt2_gt3 += gt2_gts_same_num
                    print("\t".join([b, "%s" % start, "%s" % end, "%s" % gt1_gts_same_num, "%s" % gt2_gts_same_num]))
                    #if count_diff_num == 0:
                    #    ratio =0
                    #else:
                    #    ratio =  math.log(count_diff_num,10)
                    #if ratio <= float(level[0]):
                    #    print("\t".join([b, "%s" % start, "%s" % end, "low_com","%s" % count_diff_num]))
                    #if float(level[0]) < ratio <= float(level[1]):
                    #    print("\t".join([b, "%s" % start, "%s" % end, "mid_com","%s" % count_diff_num]))
                    #if float(level[1]) < ratio :
                    #    print("\t".join([b, "%s" % start, "%s" % end, "high_com","%s" % count_diff_num]))
                # #
                start += bin_size
                end += bin_size
                gt1_gts_same_num = 0
                gt2_gts_same_num = 0

                #
            # if IsGeno(g1):
            #     count1 += 1
            # if IsGeno(g2):
            #     count2 += 1
            if (int(DP_low[0]) <= DP1 <= int(DP_high[0])
                and int(DP_low[0]) <= DP2 <= int(DP_high[0])
                and int(DP_low[0]) <= DPs <= int(DP_high[0])
                and GQ1 >= int(GQ_sample)
                and GQ2 >= int(GQ_sample)
                and GQs >= int(GQ_sample)):
                if (IsGeno(GT1) and IsGeno(GT2) and IsGeno(GTs)
                    and IsHomo(GT1) and IsHomo(GT2) and IsHomo(GTs)):
                    if GT1 != GTs and GT2 != GTs:
                        pass
                    elif GT1 == GTs and GT2 == GTs:
                        pass
                    else:
                        if GT1 == GTs:
                            gt1_gts_same_num += 1
                        if GT2 == GTs:
                            gt2_gts_same_num += 1
        #
        #deli = abs(count2-count1)
        #ratio = count3/deli
        START = str(start)
        judge = 0
        #deli = abs(count2-count1)
        for i in range(0,len(define_list)):
            if START == define_list[i][1]:
                if define_list[i][3] != 'wait_define_low':
                    #print("\t".join(str(i) for i in define_list[i][0:4]))
                #print("break1")
                    judge = 1
                else:
                    judge = 2
        if judge==2:
            count_same_gt1_gt2 += gt1_gts_same_num
            count_same_gt2_gt3 += gt2_gts_same_num
            print("\t".join([b, "%s" % start, "%s" % end, "%s" % gt1_gts_same_num, "%s" % gt2_gts_same_num]))

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
    parser.add_option("-c", dest="chr_col",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-p", dest="pos_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 2)
    parser.add_option("--DP_high", dest="DP_high",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--DP_low", dest="DP_low",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--GQ_sample", dest="GQ_sample",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--LEVEL", dest="LEVEL",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-r", dest="ref_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 3)
    parser.add_option("-a", dest="alt_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 4)
    parser.add_option("-1", dest="s1_col",
                  help="column id for 1st parent [default: %default]", metavar="INT",
                      default = 5)
    parser.add_option("-2", dest="s2_col",
                  help="column id for 2nd parent [default: %default]", metavar="INT",
                      default = 8)
    parser.add_option("-s", dest="s3_col",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 11)
    parser.add_option("--chromosome", dest="chromosome",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = '1B')
    parser.add_option("--pre_define", dest="pre_define",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
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
    GenoPhase(infile=options.infile, bin_size=int(options.bin_size), chromosome=options.chromosome,
              chr_col=int(options.chr_col), pos_col=int(options.pos_col), pre_define=options.pre_define,
              ref_col=int(options.ref_col), alt_col=int(options.alt_col), DP_high=options.DP_high, GQ_sample=options.GQ_sample,
              DP_low=options.DP_low, LEVEL=options.LEVEL,
              s1_col=int(options.s1_col), s2_col=int(options.s2_col), s3_col=int(options.s3_col))
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
