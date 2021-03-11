#!/usr/bin/env python3
# -*- coding:utf-8 -*-
def PIC ( infile="inflie.txt", del_num=0, dup_num=0, pop=100, outfile="outfile.txt"):
    with open(infile) as f:
        GSR_mat=f.readlines()
#"/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/mcl_dir/chr1A.1_mcl_undefined"
    lis=list()
    lis1=list()
    for i in GSR_mat:
        tokens=i.strip().split("\t")
        pro=np.square(len(tokens)/pop)
        lis.append(pro)
        lis1.append(len(tokens))
    rest=pop-sum(lis1)-del_num-dup_num
    rest_pro=np.square(1/pop)
    if rest!=0:
        lis=lis+[rest_pro]*rest
    if del_num != 0:
        del_pro=np.square(del_num/pop)
        lis=lis+[del_pro]
    if dup_num != 0:
        dup_pro=np.square(dup_num/pop)
        lis=lis+[dup_pro]
    print(1-sum(lis))

from optparse import OptionParser
import math
import numpy as np

# ===========================================
def main():
    parser = OptionParser()
    parser.add_option("-i", dest="infile")
    parser.add_option("-c", dest="del_num", default=0)
    parser.add_option("-d", dest="dup_num", default=0)
    parser.add_option("-p", dest="pop", default=100)
    parser.add_option("-o", dest="outfile")
    (options, args) = parser.parse_args()
    #
    PIC (options.infile,outfile=options.outfile, pop=int(options.pop), del_num=int(options.del_num),dup_num=int(options.dup_num))

# ===========================================
if __name__ == "__main__":
    main()
