#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from __future__ import division  # very important !!!!
from collections import defaultdict
import os
import gzip
import sys
import re
import math
from scipy import stats

def Count(infile, outfile, combine_data, count, compensent, normalize, split, chrom, pop_num, compensent_duplication,crossover):
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

    dic1={}
    if count == 'on' :
        for line in IN:
            line = line.strip().split('\t')
            dic1[line[1]] = dic1.setdefault(line[1], 0)+1
        for key,value in dic1.items():
            print('{key}\t{key2}\t{value}'.format(key=key, key2=int(key)+1000000, value=value/int(pop_num)), sep='\t')
#
    if crossover == 'on':
        import pandas as pd
        import numpy as np
        dic={'low':0,'mid':1,'high':1}
#        firstline = IN.readline().strip()
#        item = firstline.strip().split("\t")
#        last_pos = int(item[1])
#        last_con = dic[item[3]]
        list_1 = []
        lis = [0]
        for line in IN:
            item = line.strip("\n").split("\t")
            if "CNV" not in item[3]:
                next_pos = int(item[1])
                next_con = dic[item[3]]
                list_1.append([item[0],item[1],item[2],next_con])
        frame_1 = pd.DataFrame(list_1)
        frame_1.loc[:,1:3] = frame_1.loc[:,1:3].astype(int)
        def function_1(a,b,c,d,e,f):
            num=0
            if a[3] == b[3] == c[3] != d[3]:
                if (a[1] + b[1] + c[1]) == b[1]*3 and d[1] == c[1]+1000000:
                    num=1
            if d[3] == e[3] == f[3] != c[3]:
                if (d[1] + e[1] + f[1]) == e[1]*3 and c[1]+1000000 == d[1]:
                    num=1
            if num == 1:
                return 1
            else:
                return 0
        def function_2(a,b,c,d):
            num=0
            if a[3] == b[3] == c[3] != d[3]:
                if (a[1] + b[1] + c[1]) == b[1]*3 and abs(a[1]-d[1])==1000000:
                    num=1
            if num == 1:
                return 1
            else:
                return 0
        def function_3(a,b,c,d):
            num=0
            if a[3] == b[3] == c[3] != d[3]:
                if (a[1] + b[1] + c[1]) == b[1]*3 and abs(a[1]-d[1])==3000000:
                    num=1
            if num == 1:
                return 1
            else:
                return 0
        for i in range(1,3):
            num_1 = function_2(frame_1.loc[i,],frame_1.loc[i+1,],frame_1.loc[i+2,],frame_1.loc[i-1,])
            lis.append(num_1)
        num_2 = 0
        #for i in range(1,len(frame_1)):
        #    if frame_1.loc[i-1,][3] != frame_1.loc[i,][3] and frame_1.loc[i-1,][1] == frame_1.loc[i,][1]-1000000:
        #        lis.append(1)
        #    else:
        #        lis.append(0)
        for i in range(3,len(frame_1)-2): #len(frame_1)-4
            #lis.append(function_1(frame_1.loc[i-3,],frame_1.loc[i-2,],frame_1.loc[i-1,],frame_1.loc[i,],frame_1.loc[i+1,],frame_1.loc[i+2,],frame_1.loc[i+3,]))
            num_1 = function_1(frame_1.loc[i-3,],frame_1.loc[i-2,],frame_1.loc[i-1,],frame_1.loc[i,],frame_1.loc[i+1,],frame_1.loc[i+2,])
            if num_2 == 0 :
                lis.append(num_1)
            elif num_2 == num_2:
                lis.append(0)
            else:
                lis.append(num_1)
            num_2 = num_1
        for i in range(len(frame_1)-2,len(frame_1)):
            num_1 = function_3(frame_1.loc[i-3,],frame_1.loc[i-2,],frame_1.loc[i-1,],frame_1.loc[i,])
            lis.append(num_1)
        l4 = [item for item in lis if item == 1]
        infile_str=infile.strip("\n").split("/")
        #print(len(l4),'\t',infile_str[7])
        frame_2 = pd.DataFrame(lis)
        frame_2.columns=[4]
        df = pd.concat([frame_1,frame_2],axis=1)        
        df_new = df.loc[:,[0,1,2,4]]
        df_filter = df_new[df_new[4].isin([1])]
        df_filter.to_csv(outfile,index=False,header=None,sep='\t')

    if normalize == 'on':
        lis=[]
        for line in IN:
            item = line.strip("\n").split("\t")
            if item[0] != '0':
                lis.append(item[0])
        mode_num=float(stats.mode(lis)[0][0])
        print("%.2f" % mode_num)

    if split == 'on':
        for line in IN:
            item = line.strip("\n").split("\t")
            for i in range(int(item[1]),int(item[2]),1000000):
                print(item[0],i,i+1000000,item[3],sep="\t")
    if compensent == 'on' :
        dic1 = defaultdict(list)
        for line in IN:
            item = line.strip("\n").split("\t")
            dic1[item[0]].append(item[1])
            dic1[item[0]].append(item[2])
        wholen = (594102056,
                    689851870,
                    495453186,
                    780798557,
                    801256715,
                    651852609,
                    750843639,
                    830829764,
                    615552423,
                    744588157,
                    673617499,
                    509857067,
                    709773743,
                    713149757,
                    566080677,
                    618079260,
                    720988478,
                    473592718,
                    736706236,
                    750620385,
                    638686055)

        b = chrom
        dic2 = {"A":"0","B":"1","D":"2"}
        num = (int(b[0])-1)*3
        ch = int(dic2[b[1]])
        loc = num + ch
        if '1' in dic1:
            use_num = str(float(dic1["1"][1]))
            print("\t".join(['1', '1000001', "%s" % use_num]))
        else:
            print("\t".join(['1', '1000001', "0.0"]))
        for item in range(1, int(wholen[loc]/1000000)+1):            
            item = str(item*1000000+1)
            if item in dic1:
                use_num = float(dic1[item][1])
                print("\t".join(["%s" % item, "%s" % dic1[item][0], "%s" % use_num]))
            else:
                print("\t".join(["%s" % item, "%s" % str(int(item)+1000000), "0.0"]))


    if compensent_duplication == 'on' :
        dic1 = defaultdict(list)
        for line in IN:
            item = line.strip("\n").split("\t")
            dic1[item[0]].append(item[1])
            dic1[item[0]].append(item[2])
        wholen = (594102056,
                    689851870,
                    495453186,
                    780798557,
                    801256715,
                    651852609,
                    750843639,
                    830829764,
                    615552423,
                    744588157,
                    673617499,
                    509857067,
                    709773743,
                    713149757,
                    566080677,
                    618079260,
                    720988478,
                    473592718,
                    736706236,
                    750620385,
                    638686055)

        b = chrom
        dic2 = {"A":"0","B":"1","D":"2"}
        num = (int(b[0])-1)*3
        ch = int(dic2[b[1]])
        loc = num + ch
        if '1' in dic1:
            use_num = str(1.0+float(dic1["1"][1]))
            print("\t".join(['1', '1000001', "%s" % use_num]))
        else:
            print("\t".join(['1', '1000001', "1.0"]))
        for item in range(1, int(wholen[loc]/1000000)+2):            
            item = str(item*1000000+1)
            if item in dic1:
                use_num = 1.0+float(dic1[item][1])
                print("\t".join(["%s" % item, "%s" % dic1[item][0], "%s" % use_num]))
            else:
                print("\t".join(["%s" % item, "%s" % str(int(item)+1000000), "1.0"]))
    if infile:
        IN.close()
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
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("--combine_data", dest="combine_data",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("--count", dest="count",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--split", dest="split",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--compensent", dest="compensent",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--chrom", dest="chrom",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = '1A')
    parser.add_option("--pop_num", dest="pop_num",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = '1')
    parser.add_option("--normalize", dest="normalize",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--compensent_duplication", dest="compensent_duplication",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--crossover", dest="crossover",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
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
    Count(infile=options.infile, combine_data=options.combine_data, count=options.count, chrom=options.chrom,
          normalize=options.normalize, compensent=options.compensent, split=options.split,outfile=options.outfile,
          pop_num=options.pop_num, compensent_duplication=options.compensent_duplication,crossover=options.crossover)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
