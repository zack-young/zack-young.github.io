#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@author:yangzz
@file:comparison.py
@time:2018/9/3 17:35
"""

from __future__ import division  # very important !!!!
import os
import gzip
import sys
import re

# <a href="http://www.runoob.com/">访问菜鸟教程</a>
# <a href="https://urgi.versailles.inra.fr/wheatis/#result/term=TraesCS4B01G042200">TraesCS4B01G042200</a> 
# <a href="http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=LOC_Os03g48930">LOC_Os03g48930</a> 
# <a href="https://www.arabidopsis.org/servlets/TairObject?name=AT1G77720&type=locus">AT1G77720</a>  /data/annotation/wheat/CS_IWGSC/v1/commom/IWGSC_annotation_blastp.csv
# wheatexpress


def makehostDict(infile):
    host_dict = {}
    f_allhost = open(infile, 'r')   #rb是按二进制读取
    lines = f_allhost.readlines()
    for line in lines:
        gene = line.strip('\n')
        chrome = line[7:9]
        if gene is not '':
            if host_dict.__contains__(chrome):   # python3 不含.has_key()
                host_dict.get(chrome).append(gene) # 此处为关键向字典里已经有的key(name)值后继续添加value(host)
            else:
                host_dict.setdefault(chrome, []).append(gene)#创建{name，[host]}value为列表的格式的字典。
    #
    print('''
<table border="1">
<tr>
  <th>geneI\tDexpression\tannotation\tAth\trice</th>
</tr>'''+'\n')      
    for chrome in host_dict:
        with open("/home/wangzh/database/annotation/WheatHomo/"+ chrome, 'r') as firstdata:
            for line in firstdata.readlines()[1:]:
                line = line.strip().split('\t')
                line0 = line[0].strip('"')
                line3 = line[3].strip('"')
                line1 = line[1].strip('"')
                # 提染色体 uniq 同一染色体放一起  TraesCS5A01G42230 s[7:9]
                # theList = ['a','b','c']
                # if 'a' in theList:
                # print 'a in the list'
                if line0 in host_dict[chrome]:
                    if re.match('NA*', line[1]):
                        if re.match('NA', line[3]):
                            print('<tr>'+'\n'+'  <td><a href="https://urgi.versailles.inra.fr/wheatis/#result/term=%s">%s</a>' % (line0, line0)+ 
                                  '\t'+'<a href="http://www.wheat-expression.com/genes/show?gene_set=RefSeq1.0&name=%s&search_by=gene">%s</a>' % (line0, line0)+
                                  '\t'+
                                   + line[6] + '\t' + 'NA' + '\t'+ 'NA</td>'+ '\n' + '</tr>'+'\n')
                        else:
                            print('<tr>'+'\n'+'  <td><a href="https://urgi.versailles.inra.fr/wheatis/#result/term=%s">%s</a>' % (line0, line0)+
                                  '\t'+ '<a href="http://www.wheat-expression.com/genes/show?gene_set=RefSeq1.0&name=%s&search_by=gene">%s</a>' % (line0, line0)+
                                  '\t'+ line[6] + '\t'+'<a href="http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=%s">%s</a>' % (line1, line1)+
                                  '\t'+ 'NA</td>'+ '\n' + '</tr>'+'\n')
                    else:
                        print('<tr>'+'\n'+'  <td><a href="https://urgi.versailles.inra.fr/wheatis/#result/term=%s">%s</a>' % (line0, line0)+
                              '\t'+ '<a href="http://www.wheat-expression.com/genes/show?gene_set=RefSeq1.0&name=%s&search_by=gene">%s</a>' % (line0, line0)+
                              '\t'+ line[6] + '\t'+'<a href="http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=%s">%s</a>' % (line1, line1)+
                              '\t'+'<a href="http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=%s">%s</a></td>' % (line3, line3) + '\n' + '</tr>'+'\n')
            # else:
            #     line = line.strip().split()
                # word1 = line[0].strip().split(':')  ctrl+/
                # word2 = line[10].strip().split(':')
                # word3 = line[11].strip().split(':')
                # TODO：加制表符&emsp; 加分行符<br>
    if infile:
        f_allhost.close()


from optparse import OptionParser
# ===========================================
def main():
    usage = "Usage: %prog [-i <input>] [-o <output>]\n" \
            "Author : Yang, zhengzhao; yangzhengzhao@cau.edu.cn; 2018-09-07\n" 
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    (options, args) = parser.parse_args()
    #
    if (options.outfile is not None) :  #??????
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')
        #
    #
    makehostDict(options.infile)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
