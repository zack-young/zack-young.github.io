#!/usr/bin/env python

import os
import gzip
import sys
import re
from decimal import localcontext,Decimal

def isort (tokens, order, sample_header, debug):
    if int(debug)>0:
        sys.stderr.write("####################################\n")
    header = tokens[:3]
    samples = tokens[3:]
    order_id = []
    group = [[]]
    output_sample = ["" for i in range(0, len(sample_header))]
    #load origon sample group
    for i in range(0, len(samples)):
        if int(len(group))-1 < int(samples[i]):
            group.append([i])
        else:
            group[int(samples[i])].append(i)
    if int(debug)>0:
        sys.stderr.write("[group]\t" + str(group) + "\n")
    #load wanted sample order
    for i in range(0, len(order)):
        for j in range(0, len(sample_header)):
            if(order[i] == sample_header[j]):
                order_id.append(j)
                break
    if int(debug)>0:
        sys.stderr.write("[order_id]\t" + str(order_id) + "\n")
    #map order to sample index
    for i in range(0, len(order_id)):
        for j in range(0, len(group)):
            if order_id[i] in group[j] :
                if j == 0:
                    for k in group[j]:
                        if output_sample[k] != "":
                            break
                        output_sample[k] = 0
                    break
                else:
                    for k in group[j]:
                        if output_sample[k] != "":
                            break
                        output_sample[k] = i+1
                    break
    if int(debug)>0:
        sys.stderr.write("[output_sample]\t" + str(output_sample) + "\n")
        sys.stderr.write("####################################\n")
    #regroup
    g = len(order_id)+1
    g_add = False
    for i in range(0,len(group)):
        g_add = False
        for j in group[i]:
            if output_sample[j] != "":
                break
            if i == 0:
                output_sample[j] = 0
                continue
            output_sample[j] = g
            g_add = True
        if g_add:
            g += 1
    #output
    output_string = "%s\t%s\t%s\t" % (str(header[0]), str(header[1]), str(header[2]))
    for i in order_id:
        output_string = output_string + str(output_sample[i]) + "\t"
    for i in range(0,len(output_sample)):
        if not (i in order_id):
            output_string = output_string + str(output_sample[i]) + "\t"
    print (output_string)



def Resort ( infile, sample_order="", 
                chr_col=1, s_col=4, debug=0 ):
    #
    sys.stderr.write("Script start...\n")
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
        sys.stderr.write("\n[Error]:\n\t File cannot be open: %s\n" % infile )
        exit(-1)
    #
    # init
    sys.stderr.write("Loading data...\n")
    notcomm = False
    sample_header = []
    if sample_order == "":
        return
    order = sample_order.strip().split(",")
    for line in IN :
        if line[0:5] == "CHROM":
            notcomm = True
            sample_header = line.strip().split("\t")
            sample_header = sample_header[int(s_col)-1:]
            sample_matrix = [ [] for i in range(0,len(sample_header))]
            output_string = "%s\t%s\t%s\t" % ("CHROM", "start", "end")
            for i in order:
                if not i in sample_header:
                    sys.stderr.write("\n[Error]:\n\t Sample does not belong to input file: %s\n" % i )
                    exit(-1)
            for i in order:
                output_string = output_string + i + "\t"
            for i in sample_header:
                if not (i in order):
                    output_string = output_string + i + "\t"
            print ("%s" % output_string)
            if int(debug)>0:
                sys.stderr.write("================================\n")
                sys.stderr.write("[order]\t" + str(order) + "\n")
                sys.stderr.write("================================\n")
            continue
        #
        if notcomm:
            tokens = line.strip().split("\t") 
            isort (tokens, order, sample_header, debug)
            #
        #
    #
#

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog [-i <input file>] [-O <sample order>] [-o <output file>]\n" \
            "Author : Wang Wenxi"\
            "Update : 2018-08-11\n" \
            "Input format:\n" \
            "   VCF file\n" \
            "Output format:\n" \
            "   chr\tstart\tend\tsample1\tsample2\tsample3\n" \
            "   chr1A\t1\t1000000\t0\t0\t1\n"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-c", dest="chr_col",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-s", dest="s_col",
                  help="column id for first sample [default: %default]", metavar="INT",
                      default = 4)
    parser.add_option("-O", dest="sample_order",
                  help="Wanted order of samples", metavar="STRING",
                      default = "")
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-d", dest="debug",
                  help="weather to print more message for debuging ; "
                  , metavar="INT", default=0)
    #
    (options, args) = parser.parse_args()
    #
    if (options.outfile is not None) :
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')
        #
    #
    Resort (options.infile, options.sample_order, 
               options.chr_col, options.s_col, options.debug)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
