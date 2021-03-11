#!/data/user/yangzz/worktools/anaconda3/bin/python
###/usr/bin/env python3
# -*- coding:utf-8 -*-
def Entropy ( infile="inflie.txt", dup_num=0, del_num=0, dup_sample="a", del_sample="a",typ="file",start=1, CHR="chr1A",outfile="outfile.txt"):
    #GSR_mat=open(infile,'r') #"/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/mcl_dir/chr1A.1_mcl_undefined"
    with open(infile) as f:
        GSR_mat=f.readlines()
    lis=list()
    lis_use=list()
    lis_use_all=list()
    del_lis=del_sample.strip().split(':')
    dup_lis=dup_sample.strip().split(':')
    cnv_sample=del_lis+dup_lis
    for i in GSR_mat:
        tokens=i.strip().split("\t")
        if typ=='ALL':
            lis_use.append(tokens)
        else:
            lis_use.append([s for s in tokens if typ in s[0:3]])
    lis_use=[s for s in lis_use if len(s)>1]
    #rest=pop-sum(lis1) #-del_num-dup_num
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
        #[lis_use_all.append([i]) for i in sample_lis if typ == i[0:3]]
    sample_count=len(lis_use_all)
    lis_use=list(filter(None, lis_use))
    lis_rest1=list(set(lis_use_all).difference(set(sum(lis_use,[]))))
    lis_rest2=list(set(lis_rest1).difference(set(cnv_sample)))
    if len(lis_rest2)==0:
        lis_rest2=["none"]
    sta_hap= sample_count -del_num-dup_num-len(sum(lis_use,[]))
#--------------------------------------
    #type_sample=locals()[typ+'_sample']
    #lis_type=locals()['lis_'+typ]
    #lis2=sum(lis,[])
    #rest_lis=list(set(sample_lis).difference(set(lis2)))
    with open("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/"+CHR+"_"+typ+"_hap_count.txt", 'a+') as f:
        hap_num=1
        for i in lis_use:
            if len(i)!=1:
                f.write('\t'.join([CHR,str(start),str(len(i)),str(hap_num),','.join(i)]))
                hap_num+=1
                f.write("\n")
        f.write('\t'.join([CHR,str(start),str(sta_hap),"exclusive",','.join(lis_rest2)]))
        f.write("\n")
        f.write('\t'.join([CHR,str(start),str(del_num),"deletion",','.join(del_lis)]))
        f.write("\n")
        f.write('\t'.join([CHR,str(start),str(dup_num),"duplication",','.join(dup_lis)]))
        f.write("\n")
        
#    with open("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/landrace_hap.txt", 'a+') as f:
#        f.write('\n'.join([str(j) for j in [len(s) for s in lis_L]]))
#        f.write("\n")    
    hap_count_2=0
    hap_count_1=0
    hap_count_rest=0
    C_hap_count_2=0
    C_hap_count_1=0
    C_hap_count_rest=0
    L_hap_count_2=0
    L_hap_count_1=0
    L_hap_count_rest=0
    sample_count_2=0
    #for a in range(0,len(lis)):
    #    i=len(lis[a])
    #    if i==1:
    #        hap_count_1+=1
    #    if i>=2 and i/all_hap < 0.05:
    #        hap_count_2+=1
    #    if i/all_hap >= 0.05:
    #        hap_count_rest+=1
    #for a in range(0,len(lis_C)):
    #    i=len(lis_C[a])
    #    if i==1:
    #        C_hap_count_1+=1
    #    if i>=2 and i/C_hap < 0.05:
    #        C_hap_count_2+=1
    #    if i/C_hap >= 0.05:
    #        C_hap_count_rest+=1
    #for a in range(0,len(lis_L)):
    #    i=len(lis_L[a])
    #    if i==1:
    #        L_hap_count_1+=1
    #    if i>=2 and i/L_hap < 0.05:
    #        L_hap_count_2+=1
    #    if i/L_hap >= 0.05:
    #        L_hap_count_rest+=1
    #with open("/data/user/yangzz/mapping/fieldergenomecompare/statistic/GSR_MCL/"+CHR+"1_2_CN_hap.txt", 'a+') as f:
    #    f.write('\t'.join([CHR, "%s" % start, "%s" % hap_count_1 , "%s" % hap_count_2,"%s" % hap_count_rest,"%s" % C_hap_count_1 , "%s" % C_hap_count_2, "%s" % C_hap_count_rest,"%s" % L_hap_count_1 , "%s" % L_hap_count_2,"%s" % hap_count_rest]))
    #    f.write("\n")
#    if del_num != 0:
#        del_pro=del_num/pop
#        del_shannonEnt=-del_pro * np.math.log(del_pro, 2)
#        lis=lis+[del_shannonEnt]
#    if dup_num != 0:
#        dup_pro=dup_num/pop
#        dup_shannonEnt=-dup_pro * np.math.log(dup_pro, 2)
#        lis=lis+[dup_shannonEnt]
    #print(sum(lis))

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
