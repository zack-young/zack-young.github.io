#!/data/user/yangzz/worktools/anaconda3/bin/python
###/usr/bin/env python3
# -*- coding:utf-8 -*-
def Entropy ( infile="inflie.txt", dup_num=0, del_num=0, dup_sample="a", del_sample="a",typ="file",start=1, CHR="chr1A",outfile="outfile.txt"):
    #GSR_mat=open(infile,'r') #"/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/mcl_dir/chr1A.1_mcl_undefined"
    sample_lis=[]
    with open(infile) as f:
        GSR_mat=f.readlines()
    with open("/data/user/yangzz/mapping/fieldergenomecompare/202009_11_yaoyy/array_quality/good_q_sample.txt") as f:
        sample=f.readlines()
    for i in sample:
        tokens=i.strip().split("\t")
        sample_lis.append(tokens[0])
    lis=list()
    lis_use=list()
    for i in GSR_mat:
        tokens=i.strip().split("\t")
        sam=tokens[4].split(',')
        if tokens[3]=='exclusive':
            res_sam=list(set(sample_lis).intersection(set(sam)))
            for item in res_sam:
                print('\t'.join(['\t'.join([str(j) for j in tokens[0:3]]),'exclusive',item,'1']))
        else:
            print('\t'.join(tokens))
        
from optparse import OptionParser
import math
import numpy as np

# ===========================================
def main():
    parser = OptionParser()
    parser.add_option("-i", dest="infile")
    parser.add_option("-c", dest="del_num", default=0)
    parser.add_option("-d", dest="dup_num", default=0)
    parser.add_option("--del_sample", dest="del_sample", default=0)
    parser.add_option("--dup_sample", dest="dup_sample", default=0)
    parser.add_option("-t", dest="typ", default=100)
    parser.add_option("-s", dest="start", default=100)
    parser.add_option("-r", dest="CHR", default=100)
    parser.add_option("-o", dest="outfile")
    (options, args) = parser.parse_args()
    #
    Entropy (infile=options.infile, del_num=int(options.del_num),dup_num=int(options.dup_num),del_sample=options.del_sample,dup_sample=options.dup_sample,typ=options.typ,start=int(options.start),CHR=options.CHR , outfile=options.outfile)

# ===========================================
if __name__ == "__main__":
    main()
