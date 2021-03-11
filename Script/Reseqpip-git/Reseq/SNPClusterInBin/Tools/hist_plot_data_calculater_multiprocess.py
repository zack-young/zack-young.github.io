#!/usr/bin/env python

import os
import gzip
import sys
import re
import multiprocessing
from decimal import localcontext,Decimal

def calculate(bin_size, s1, s2, debug):
    count = 0
    for i in range(0,len(s1)):
        if s1[i] != s2[i]:
            count += 1
        #
    #
    dec_count = Decimal(str(count))
    dec_bin_size = Decimal(str(bin_size))
    if int(debug) > 2:
        print ("[cal_s1, len(s1)]\t%s,\t%d", [str(s1), len(s1)])
        print ("[cal_s2, len(s2)]\t%s,\t%d", [str(s2), len(s2)])
        print ("[count, value]\t%d,\t%f", [str(count), str(dec_count/dec_bin_size)])
    #
    return dec_count/dec_bin_size
#

def analyze(chr_cur, start, end, bin_size, sample_matrix, debug, writing_lock):
    if int(debug) > 1:
        print("[*] call analyze")
    #
    if(len(sample_matrix[0]) == 0):
        return
    #
    loop_swi = True
    for i in range(0, len(sample_matrix)-1):
        for j in range(i+1, len(sample_matrix)):
            v = calculate(bin_size, sample_matrix[i], sample_matrix[j], debug)
            while loop_swi:
                if not writing_lock[0]:
                    writing_lock[0] = True
                    loop_swi = False
            print(v)
            writing_lock[0] = False
        #
    #
#

def call_analyze(chr_cur, start, end, bin_size, sample_matrix, debug, process_status, process_id, writing_lock):
    process_status[process_id] = False
    analyze(chr_cur, start, end, bin_size, sample_matrix, debug, writing_lock)
    process_status[process_id] = True
#

def Cluster ( infile, bin_size = 1000000, max_pros = 30,
                chr_col=1, pos_col=2, s_col=10, debug=0 ):
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
    start = 1
    end = start + bin_size - 1
    notcomm = False
    chr_cur = ""
    sample_header = []
    sample_matrix = []

    process_status = [True for i in range(0,max_pros)]
    writing_lock = [False]
    loop_swi = True

    for line in IN :
        if int(debug) > 1:
            print ("========================================")
            print ("[line]\t%s" % str(line[:20])) 
            print ("[chr_cur, start, end, notcomm, debug]\t%s\t%s\t%s\t%s\t%s" % [str(chr_cur), str(start), str(end), str(notcomm), str(debug)])
        #
        if int(debug) > 2:
            print (str(sample_header))
            print (str(sample_matrix))
        #
        if line[0:6] == "#CHROM":
            notcomm = True
            sample_header = line.strip().split("\t")
            sample_header = sample_header[int(s_col)-1:]
            sample_matrix = [ [] for i in range(0,len(sample_header))]
            continue
        #
        if notcomm:
            tokens = line.strip().split("\t") 
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            if int(debug) > 1:
                print ("[chr, pos]\t%s\t%s" % [str(chr), str(pos)])
                print ("========================================")
            if chr == chr_cur :
                while pos > end :
                    if int(debug) > 1:
                        print ("========================================")
                        print ("[loop]", "001")
                        print ("========================================")
                    #thread start
                    while loop_swi:
                        for i in range(0, len(process_status)):
                            if process_status[i]:
                                pro = multiprocessing.Process(target=call_analyze, args=[chr_cur, start, end, bin_size, sample_matrix, debug, process_status, i, writing_lock])
                                pro.start()
                                loop_swi = False
                                break
                            #
                        #
                    #
                    loop_swi = True
                    start += bin_size
                    end += bin_size
                    sample_matrix = [ [] for i in range(0,len(sample_header))]
                #
            else : # Start a new chromosome
                if chr_cur is not "" :
                    #thread start
                    while loop_swi:
                        for i in range(0, len(process_status)):
                            if process_status[i]:
                                pro = multiprocessing.Process(target=call_analyze, args=[chr_cur, start, end, bin_size, sample_matrix, debug, process_status, i, writing_lock])
                                pro.start()
                                loop_swi = False
                                break
                            #
                        #
                    #
                    loop_swi = True
                #
                chr_cur = chr
                start = 1
                end = start + bin_size - 1
                sample_matrix = [ [] for i in range(0,len(sample_header))]
                while pos > end :
                    if int(debug) > 2:
                        print ("========================================")
                        print ("[loop, pos, end]\t%s\t%s\t%s" % ["002", str(pos), str(end)])
                        print ("========================================")
                    #
                    start += bin_size
                    end += bin_size
                    sample_matrix = [ [] for i in range(0,len(sample_header))]
                    #thread start
                    while loop_swi:
                        for i in range(0, len(process_status)):
                            if process_status[i]:
                                pro = multiprocessing.Process(target=call_analyze, args=[chr_cur, start, end, bin_size, sample_matrix, debug, process_status, i, writing_lock])
                                pro.start()
                                loop_swi = False
                                break
                            #
                        #
                    #
                    loop_swi = True
                #
            #
            line_elements = line.strip().split("\t")
            for i in range(0,len(sample_matrix)):
                if (line_elements[i + s_col -1])[0:3] == "1/1":
                    sample_matrix[i].append(1)
                else:
                    sample_matrix[i].append(0)
                #
            #
        #
    #
    #thread start
    while loop_swi:
        for i in range(0, len(process_status)):
            if process_status[i]:
                pro = multiprocessing.Process(target=call_analyze, args=[chr_cur, start, end, bin_size, sample_matrix, debug, process_status, i, writing_lock])
                pro.start()
                loop_swi = False
                break
            #
        #
    #
    loop_swi = True
    while loop_swi:
        loop_swi = False
        for i in process_status:
            if not i:
                loop_swi = True
                break
            #
        #
    #
    if infile :
        IN.close()
    #
#

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog [-i <input file>] [-b 1000000] [-C 0.001] [-o <output file>]\n" \
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
    parser.add_option("-b", dest="bin_size",
                  help="size of bin (bp) [default: %default]", metavar="INT",
                      default = 1000000)
    parser.add_option("-c", dest="chr_col",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-P", dest="max_pros",
                  help="maxium of processes [default: %default]", metavar="INT",
                      default = 30)
    parser.add_option("-p", dest="pos_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 2)
    parser.add_option("-s", dest="s_col",
                  help="column id for first sample [default: %default]", metavar="INT",
                      default = 10)
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
    Cluster (options.infile, int(options.bin_size), int(options.max_pros), 
               options.chr_col, options.pos_col, options.s_col, options.debug)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
