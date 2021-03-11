#!/usr/bin/env python
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  C4-2    JM22t   M1-33   M1-37   M61-1   S138    Y968    s3753-1
#chr1A   349     .       G       A       318.91  PASS    AC=4;AF=0.25;AN=16;DP=71;ExcessHet=0.0335;FS=0;MLEAC=4;MLEAF=0.25;MQ=59.34;QD=31.89;SOR=1.085;ANN=A|intergenic_region|MODIFIER|CHR_START-TraesCS1A01G000100|CHR_START-TraesCS1A01G000100|intergenic_region|CHR_START-TraesCS1A01G000100|||n.349G>A||||||        GT:AD:DP:GQ:PL  0/0:7,0:7:21:0,21,260   0/0:12,0:12:33:0,33,495 0/0:16,0:16:48:0,48,625 0/0:12,0:12:33:0,33,495 0/0:9,0:9:27:0,27,338   0/0:5,0:5:12:0,12,180   1/1:0,6:6:18:239,18,0   1/1:0,4:4:12:135,12,0
#17 lie
#!/usr/bin/env python
import os
import gzip
import sys
#import re
import numpy as np

def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)

def callTotalVariation(infile,outfile,wt=10,mut=11,qvalue = 30,mindepth=3):
	data = open(infile)
	out = open(outfile,'w')
	da = data.readlines()
	wt = int(wt)
	mut = int(mut)
	mind = int(mindepth)
	qvalue = float(qvalue)
	for i in da:
        	if not i.startswith('#'):
                	i =i.split('\t')
                	q = float(i[5].strip())
                	MUT = str(i[mut].strip())
               		WT = str(i[wt].strip())
                	MUT = MUT.split(':')
                	WT = WT.split(':')
                	gtY = str(MUT[0].strip())
                	gtN = str(WT[0].strip())
                	if gtY != './.' and gtN != './.':
                        	totaldepY = int(MUT[2].strip())
                        	totaldepN = int(WT[2].strip())
                        	vdepY = MUT[1].strip().split(',')
                        	ref1 = int(vdepY[0].strip())
                        	alt1 = int(vdepY[1].strip())
                        	vdepN = WT[1].strip().split(',')
                        	ref2 = int(vdepN[0].strip())
                        	alt2 = int(vdepN[1].strip())
                        	ref = i[3].strip()
                        	alt = i[4].strip()
                        	#if q > qvalue and ((gtY == '0/0' and gtN == '1/1') or (gtY == '1/1' and gtN == '0/0')) and 100 > totaldepY > mind and 100 > totaldepN > mind :
                                chr1 = i[0].strip()
                                pos = i[1].strip()
                                #ref = i[3].strip()
                                #alt = i[4].strip()
                                info = i[7].strip()
                                out.write(chr1+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+str(ref2)+'\t'+str(alt2)+'\t'+str(ref1)+'\t'+str(alt1)+'\t'+info+'\n')

from optparse import OptionParser

def main():
	usage = 'Usage: python %prog [-i <input>] [-w 10] [-m 11] [-q 30] [-d 3] [-o <output>] \n'\
		'Author: Chen,Yongming; chen_ym@cau.edu.cn; 2018-12-07 \n'\
		'Description: Get the homozygous mutations between wildtype and mutant in vcf(Wang ZH). \n'\
		'By default,input should be vcf format.\n'
	parser = OptionParser(usage)
	parser.add_option('-i',dest = 'infile',
			help='Input file' ,metavar = 'VCF')
	parser.add_option('-w',dest = 'wt',
			help='The column of wildtype in inputfile [default: %default]',default = 10,metavar = 'INT')
	parser.add_option('-m',dest = 'mut',
			help='The column of mutant in outputfile[default: %default]',default = 11,metavar = 'INT')
	parser.add_option('-o',dest = 'outfile',
			help='Output file',metavar = 'TXT')
	parser.add_option('-q',dest = 'var',
			help='The smallest quantity value of variation[default: %default]',default = 0,metavar = 'FLOAT')
	parser.add_option('-d',dest = 'dep',
			help='the smallest reads depth[default: %default]',default = 0,metavar = 'INT')
	(options,args) = parser.parse_args()
	callTotalVariation(options.infile,options.outfile,int(options.wt)-1,int(options.mut)-1,float(options.var),int(options.dep))
# ===========================================
if __name__ == "__main__":
    main()
