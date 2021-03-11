#!/data/user/yangzz/worktools/anaconda3/bin/python
## #!/usr/bin/env python3
# -*- coding:utf-8 -*-
def Entropy ( infile="inflie.txt", dup_num=0, del_num=0, pop=100, typ='all',outfile="outfile.txt"):
    with open(infile) as f:
        GSR_mat=f.readlines()
    lis=list()
    lis_use=list()
    lis_use_all=list()
    sample_lis=[]
    with open("/data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_ALL.txt") as f:
        sample=f.readlines()
    for i in sample:
        tokens=i.strip().split("\t")
        sample_lis.append(tokens[3]+"_"+tokens[2])
    if typ=='ALL':
        lis_use_all=sample_lis
    else:
        lis_use_all=[i for i in sample_lis if typ in i[0:3]]

    for i in GSR_mat:
        tokens=i.strip().split("\t")
        if typ=='ALL':
            lis_use.append(tokens)
        else:
            lis_use.append([s for s in tokens if typ in s[0:3]])

    lis_use=[s for s in lis_use if len(s)>1]
    sample_count=len(lis_use_all)
    lis_use=list(filter(None, lis_use))
    sta_hap= sample_count -del_num-dup_num-len(sum(lis_use,[]))
    use_sample = sample_count -del_num-dup_num
    #----------------------
    if use_sample != 0:
        for i in lis_use:
            pro=len(i)/use_sample
            shannonEnt=-pro * np.math.log(pro, 2)
            lis.append(shannonEnt)
        rest_pro=1/use_sample
        rest_shannonEnt=-rest_pro * np.math.log(rest_pro, 2)
        if sta_hap!=0:
            lis=lis+[rest_shannonEnt]*sta_hap
    else:
        lis=[0]
    #if del_num != 0:
    #    del_pro=del_num/pop
    #    del_shannonEnt=-del_pro * np.math.log(del_pro, 2)
    #    lis=lis+[del_shannonEnt]
    #if dup_num != 0:
    #    dup_pro=dup_num/pop
    #    dup_shannonEnt=-dup_pro * np.math.log(dup_pro, 2)
    #    lis=lis+[dup_shannonEnt]
    print(sum(lis))

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
    parser.add_option("-t", dest="typ", default='all')
    parser.add_option("-o", dest="outfile")
    (options, args) = parser.parse_args()
    #
    Entropy (options.infile, del_num=int(options.del_num), dup_num=int(options.dup_num),pop=int(options.pop), typ=options.typ, outfile=options.outfile)

# ===========================================
if __name__ == "__main__":
    main()
