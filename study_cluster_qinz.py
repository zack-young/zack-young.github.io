#!/bin/env python
# coding: utf-8
# qinz 2019.4.30 qz19940718@gmail.com
# "make ensemble-biomart out homologs gene to cluster" 

""" AET0Gv20000400	TraesCS1B02G186000
AET0Gv20001400	
AET0Gv20002000	TraesCS4A02G058800
AET0Gv20002000	TraesCS4B02G245900
AET0Gv20002300	TraesCS2A02G244100
AET0Gv20002300	TraesCS2B02G263000
AET0Gv20003000	TraesCS2A02G240300
AET0Gv20003000	TraesCS2B02G258700
AET0Gv20003100	TraesCS2A02G237900
AET0Gv20003100	TraesCS2B02G257400
AET0Gv20003200	
AET0Gv20003400	TraesCS3B02G533800
AET0Gv20003400	TraesCS3B02G533900
AET0Gv20003400	TraesCSU02G051900  """

import re
from sys import argv
outdir = {}
filelist = []
filelist2 = []
with open (argv[1]) as infile :
    for i in infile :
        i2 = i.strip().split("\t")
        if len(i2) > 1 :
            if not ((i2[0] in filelist) and (i2[1] in filelist)) :
                filelist.append(i2[0])
                filelist.append(i2[1])
        elif len(i2) == 1 :
            filelist2.append(i2[0])
            filelist2.append("")

def mergelist(list1,list2):
    for i in list2:
        if i not in list1:
            list1.append(i)
        else :
            continue
    return list1
locate = lambda x : [i for i,v in enumerate(filelist) if v==x]
def multilocate(poslist_in) :
    newNeedLocatePos = []
    NeedLocatePos = poslist_in
    for i in NeedLocatePos :
        newNeedLocatePos = mergelist(newNeedLocatePos,locate(filelist[i]))
    return newNeedLocatePos

def getcluster():

    cluster = []
    for l in range(0,len(filelist),2):
        i = filelist[l]
        poslist = [[],[]]
        Lpos = locate(i)
        Rpos = [x+1 for x in Lpos ]
        poslist[0] = mergelist(poslist[0],Lpos)
        poslist[1] = mergelist(poslist[1],Rpos)
        poslist[1] = multilocate(poslist[1])
        while len(poslist[0]) != len(poslist[1]) :
            if len(poslist[0]) < len(poslist[1]):
                poslist[0] = multilocate([i-1 for i in poslist[1]])

            elif len(poslist[1]) < len(poslist[0]):
                poslist[1] = multilocate([i+1 for i in poslist[0]])

            elif len(poslist[0]) == len(poslist[1]) :

                break
        if mergelist(poslist[0],poslist[1]) not in cluster :
            cluster.append(mergelist(poslist[0],poslist[1]))
    return cluster
clusters = getcluster()

def uniqlist(oldlist):
    newlist = list(set(oldlist))
    newlist.sort(key=oldlist.index)
    return newlist
def outcluster(clusterlist_in,setlabelist = ["TraesCS\dA\d+\w\d+","TraesCS\dB\d+\w\d+","AET\d\w+\d+","TraesCSU\d+\w\d+"]) :

    c = 0
    for i in clusterlist_in :
        c += 1
        outlist = ["cluster%s" %c,"","","","",] 
        percluster = [filelist[x] for x in i]
        for x in range(1,5) :
            outlist[x] = re.findall(setlabelist[x-1],str(percluster))
            outlist[x] = uniqlist(outlist[x])
            ratio = [len(outlist[1]),len(outlist[2]),len(outlist[3]),len(outlist[4])]
        print outlist[0] +"\t"+ ";".join(outlist[1])+"\t"+";".join(outlist[2])+"\t"+";".join(outlist[3])\
            +"\t"+";".join(outlist[4])+"\t"+":".join([str(i) for i in ratio]) +"\t"+str(sum(ratio))
    for i in filelist2 :

        c += 1
        noclulist = ["nocluster%s" %c,"","","","",]
        for x in range(1,5) :
            noclulist[x] = re.findall(setlabelist[x-1],str(i))
            noclulist[x] = uniqlist(noclulist[x])
            ratio = [len(noclulist[1]),len(noclulist[2]),len(noclulist[3]),len(noclulist[4])]
        if sum(ratio) != 0 :
            print noclulist[0] +"\t"+ ";".join(noclulist[1])+"\t"+";".join(noclulist[2])+"\t"+";".join(noclulist[3])\
            +"\t"+";".join(noclulist[4])+"\t"+":".join([str(i) for i in ratio]) +"\t"+str(sum(ratio)) 

outcluster(clusters)
















