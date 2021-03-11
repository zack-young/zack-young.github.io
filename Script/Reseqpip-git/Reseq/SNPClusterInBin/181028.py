#!/usr/bin/env python

import os
import gzip
import sys
import re
from copy import deepcopy
from decimal import localcontext,Decimal

def main_sample_sel (header, data, m_b):
    main_sam_lis = []
    black_block = (len(header)-3)*(len(data))
    #black_percent = 0.1
    sys.stderr.write("Finding main samples ......\n")

    while(not black_block < m_b):
        max_new_sam_id = 0
        for i in range(3, len(header)):
            temp_sam_id = i
            temp_black_block = 0
            if temp_sam_id in main_sam_lis:
                continue
            for j in data:
                for k in range(3, len(j)):
                    add_switch = True
                    for x in main_sam_lis:
                        if j[k] == j[x]:
                            add_switch = False
                    if j[k] == j[temp_sam_id]:
                        add_switch = False
                    if add_switch:
                        temp_black_block += 1
            #sys.stderr.write(str(max_new_sam_id))
            #sys.stderr.write(str(temp_black_block))
            if temp_black_block < black_block:
                black_block = temp_black_block
                max_new_sam_id = temp_sam_id
        if max_new_sam_id != 0:
            main_sam_lis.append(max_new_sam_id)
        sys.stderr.write(str(main_sam_lis)+"\n")

    sys.stderr.write("Finding main samples ......Done\n")
    return main_sam_lis

def reorder_main_sam (header, data, main_sample_list):
    sys.stderr.write("Reorder main samples ......\n")
    loop_switch = True
    N1_l = deepcopy(main_sample_list)
    N1 = count_N(deepcopy(data), N1_l)
    while loop_switch:
        loop_switch = False
        for i in range(0, len(main_sample_list)-1):
            N2_l = deepcopy(N1_l)
            N2_l[i], N2_l[i+1] = N2_l[i+1], N2_l[i]
            N2 = count_N(deepcopy(data), N2_l)

            if N2 < N1:
                N1 = N2
                N1_l = deepcopy(N2_l)
                loop_switch = True
        sys.stderr.write(str(N1_l) + ", " + str(N1) + " jumps\n")
    N2 = count_N(data, N1_l)
    sys.stderr.write("Reorder main samples ......Done\n")
    return N1, N1_l

def count_N (data, main_sample_list):
    jump_count = 0
    # Deal with other blocks
    for i in data:
        for j in range(3, len(i)):
            black_switch = True
            if (j in main_sample_list) or (i[j]=="0"):
                continue
            for k in range(0, len(main_sample_list)):
                if i[j] == i[main_sample_list[k]]:
                    i[j] = -(k+1)
                    black_switch = False
            if black_switch:
                i[j] = 1
    # deal with main samples
    for j in range(0, len(main_sample_list)):
        for i in data:
            if int(i[main_sample_list[j]]) > 0:
                for k in range(j+1, len(main_sample_list)):
                    if int(i[main_sample_list[k]]) > 0:
                        if i[main_sample_list[k]] == i[main_sample_list[j]]:
                            i[main_sample_list[k]] = -(j+1)
                i[main_sample_list[j]] = -(j+1)
    # count jump number
    for i in range(3, len(data[0])):
        for j in range(0, len(data)-1):
            if data[j][i] != data[j+1][i]:
                jump_count += 1

    return jump_count


def data_read ( infile, chr_col=1, s_col=4, m_b=150 ):
    #
    sys.stderr.write("Script start......\n")
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
    sys.stderr.write("Loading data......\n")
    notcomm = False
    sample_header = []
    data_matrix = []
    for line in IN :
        if line[0:5] == "CHROM":
            notcomm = True
            sample_header = line.strip().split("\t")
            print "\t".join(sample_header)
            continue
        #
        if notcomm:
            tokens = line.strip().split("\t") 
            data_matrix.append(tokens)
            #
        #
    sys.stderr.write("Loading data......Done\n")
    main_sample_list = main_sample_sel(sample_header, data_matrix, m_b)
    jump_count, reordered_main_sample_list = reorder_main_sam(sample_header, data_matrix, main_sample_list)
    for line in data_matrix:
        print "\t".join([str(i) for i in line])
        print "\n"
    sys.stderr.write("Main sample list(reordered):\n")
    sys.stderr.write(",".join([sample_header[i] for i in reordered_main_sample_list]) + "\n")
    #
#

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog [-i <input file>] [-o <output file>] [-s <column number of first sample>]\n" \
            "Author : Wang Wenxi"\
            "Update : 2018-08-11\n" \
            "Input format:\n" \
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
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-b", dest="max_black",
                  help="max quantity of black blocks ; "
                  , metavar="INT", default=150)
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
    data_read(options.infile, options.chr_col, options.s_col, options.max_black)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
