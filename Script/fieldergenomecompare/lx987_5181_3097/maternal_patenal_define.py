#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import re
import pandas as pd
def Count(infile_1, infile_2, outfile, sample_1, sample_2, son_define, son_norm):
    # open the infile
    #
    if son_define == "on":
        with open(infile_1, 'r') as IN:
            list_1 = []
            for line in IN: 
                tokens = line.strip().split("\t")
                list_1.append(tokens)
        with open(infile_2, 'r') as IN:
            list_2 = []
            for line in IN: 
                tokens = line.strip().split("\t")
                list_2.append(tokens)
        frame_1 = pd.DataFrame(list_1)
        frame_2 = pd.DataFrame(list_2)
        df = pd.concat( [frame_1, frame_2.loc[:,3]], axis=1 )
        df.columns=['chr','start','end',sample_1,sample_2]
        df_snp = df[~(df[sample_1].str.contains('CNV'))&~(df[sample_2].str.contains('CNV'))].copy()
        df_CNV = df[(df[sample_1].str.contains('CNV'))|(df[sample_2].str.contains('CNV'))].copy()
        df_snp['belong']=df_snp[[sample_1,sample_2]].apply(lambda x: "both" if (x[sample_1] == x[sample_2] and x[sample_1] == 'low') else sample_1 if x[sample_1] == 'low'
                                                           else sample_2 if x[sample_2] == 'low' else "unknown",axis=1)
        df_CNV['belong']=df_CNV[[sample_1,sample_2]].apply(lambda x: "both" if x[sample_1] == x[sample_2]  else sample_1 if x[sample_1] == 'low' else sample_2 if x[sample_2] == 'low' 
                                                           else sample_1 if x[sample_1].find("both_CNV")>=0
                                                           else sample_2 if x[sample_2].find("both_CNV")>=0  
                                                           else "unknown" ,axis=1)
        df_new=df_snp.append(df_CNV)
    
        df_new.to_csv(outfile,index=False,sep='\t')
    if son_norm == "on":
        with open(infile_1, 'r') as IN:
            list_1 = []
            for line in IN: 
                tokens = line.strip().split("\t")
                list_1.append(tokens)
        frame_1 = pd.DataFrame(list_1)
        lis=frame_1[frame_1[3]=="both"].index.tolist()
        for i in lis:
            i=int(i)
            if i == 0:
                if frame_1.iloc[i+1,3] != "unknown":
                    frame_1.iloc[i,3]=frame_1.iloc[i+1,3]
            elif i == len(frame_1)-1:
                if frame_1.iloc[i-1,3] != "unknown":
                    frame_1.iloc[i,3]=frame_1.iloc[i-1,3]
            elif frame_1.iloc[i-1,3] == frame_1.iloc[i+1,3]:
                frame_1.iloc[i,3]=frame_1.iloc[i-1,3]    
        frame_1.columns=['chr','start','end','belong']
        frame_1.to_csv(outfile,index=False,sep='\t')




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
    parser.add_option("--i1", dest="infile_1",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("--i2", dest="infile_2",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")

    parser.add_option("--s1", dest="sample_1",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'lx987')
    parser.add_option("--son_define", dest="son_define",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--son_norm", dest="son_norm",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--s2", dest="sample_2",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = '3097')    
    #
    (options, args) = parser.parse_args()
    #
        #
    #
    Count(infile_1=options.infile_1, infile_2=options.infile_2, outfile=options.outfile, son_norm=options.son_norm,son_define=options.son_define, sample_1=options.sample_1, sample_2=options.sample_2)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
