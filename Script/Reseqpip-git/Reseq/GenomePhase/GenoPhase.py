#!/usr/bin/env python

import os
import gzip
import sys
import re

def error(msg):
    print(sys.stderr, 'ERROR: %s' % msg)
    exit(1)
#
def IsGeno ( geno ) :
    if len(geno)==3 and geno[1]=='/' and geno!="./." :
        return 1
    else :
        return 0
    #
#
def IsHete ( geno ) :
    if geno[0] == geno[2] :
        return 0
    else :
        return 1
    #
#

def IsValid (g1, g2, gs) :
    if g1 == g2 : #and g2 == gs :
        return 0
    #
    if IsGeno(g1) and IsGeno(g2) and IsGeno(gs) :
        if IsHete(g1) or IsHete(g2) or IsHete(gs) :
            return 0
        else :
            return 1
        #
    else :
        return 0
    #
#

def OutputLinehomo(chr_cur, start, end, count, count1, count2) :
    if count1 == count2 :
        state = "NA"
    elif count1 > count2 :
        state = "1"
    else :
        state = "2"
    #
    print("\t".join([chr_cur, "%s" % start, "%s" % end,
                     "%s" % count, "%s" % count1, "%s" % count2, state]) )
    #
#


def GenoPhase ( infile, bin_size,
                chr_col=1, pos_col=2, p1_col=3, p2_col=4, s_col=5 ):
    # open the infile
    try:
        if infile :
            if infile.endswith(".gz") :
                IN = gzip.open(infile, 'rb')
            else :
                IN = open(infile, 'r')
            #
        else :
            IN = sys.stdin
        #
    except IOError:
        print("\n[Error]:\n\t File cannot be open: %s" % infile)
        exit(-1)
    #
    # init
    chr_cur = None
    start = 1
    end = start + bin_size - 1
    count = 0
    count1 = 0
    count2 = 0
    #
    for line in IN :
        tokens = line.strip().split("\t")
        chr = tokens[chr_col - 1]
        pos = int(tokens[pos_col - 1])
        g1 = tokens[p1_col - 1]
        g2 = tokens[p2_col - 1]
        gs = tokens[s_col - 1]
        #
        if IsValid(g1, g2, gs) :
            if chr == chr_cur :
                while pos > end :
                    OutputLinehomo(chr_cur, start, end, count, count1, count2)
                    start += bin_size
                    end += bin_size
                    count = 0
                    count1 = 0
                    count2 = 0
                #
            else : # Start a new chromosome
                if chr_cur is not None :
                    OutputLinehomo(chr_cur, start, end, count, count1, count2)
                #
                chr_cur = chr # Beginning
                start = 1
                end = start + bin_size - 1
                while pos > end :
                    start += bin_size
                    end += bin_size
                    OutputLinehomo(chr_cur, start, end, 0, 0, 0)
                #
                count = 0
                count1 = 0
                count2 = 0
                #
            #
            # compare genotype
            count += 1
            if g1 == gs :
                count1 += 1
            else :
                count2 += 1
            #
        #
        #
    #
    # output last bin info
    OutputLinehomo(chr_cur, start, end, count, count1, count2)
    #
    if infile :
        IN.close()
    #
#



from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog [-i <input>] [-b 1000] [-o <output>]\n" \
            "Author : Guo, Weilong; guoweilong@126.com; 2018-06-07\n" \
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
    parser.add_option("-1", dest="p1_col",
                  help="column id for 1st parent [default: %default]", metavar="INT",
                      default = 5)
    parser.add_option("-2", dest="p2_col",
                  help="column id for 2nd parent [default: %default]", metavar="INT",
                      default = 7)
    parser.add_option("-s", dest="s_col",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 6)
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    #
    (options, args) = parser.parse_args()
    #
    if (options.outfile is not None) :  #??????
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')
        #
    #
    GenoPhase (options.infile, int(options.bin_size),
               options.chr_col, options.pos_col,
               options.p1_col, options.p2_col, options.s_col)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
