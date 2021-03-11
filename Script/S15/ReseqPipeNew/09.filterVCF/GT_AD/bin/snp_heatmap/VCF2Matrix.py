#!/usr/bin/env python

from optparse import OptionParser
import csv
import sys
import gzip
import copy

def VCF2Matrix (IN, codelist = "-1,0,0.5,1") :
    #
    [code_missing, code_nosnp, code_hete, code_homo] = codelist.split(",")
    #
    for line in IN :
        if line[0] != '#' :
            tokens = line.strip().split("\t")
            [CHR, POS, _, REF, ALT] = tokens[0:5]

            # in some row that has * in ALT col, snpEFF won't annotation that
            try :
                [A2,A3,A8,A11] = list( tokens[7].split("ANN=")[1].split("|")[i] for i in [1,2,7,10] )
                Tojoin = (POS, A2, A3, A8, A11)
                outline = "|".join(Tojoin)
                Ntoken = len(tokens)
                for i in range(9, Ntoken) :
                    geno = tokens[i].split(":")[0]
                    if geno == "./." :
                        outline = outline + "\t" + code_missing # NOT SNP
                    elif geno == "0/0":
                        outline = outline + "\t" + code_nosnp # NOT SNP
                    elif geno[0] == geno[2] :
                        outline = outline + "\t" + code_homo # homozygous SNP
                    else :
                        outline = outline + "\t" + code_hete # heterozygous SNP
                    #
                #
                sys.stdout.write(outline+"\n")
            except IndexError :
                pass
        else :
            if line[1] != '#' :
                tokens = line.strip().split("\t")
                Ntoken = len(tokens)
                outline  = "ANN"
                for i in range(9, Ntoken) :
                    outline = outline + "\t" + tokens[i]
                #
                sys.stdout.write(outline+"\n")
            #
        #
    #
    if IN is not sys.stdin:
        IN.close()
#   #
#


def main():
    usage = "Usage: SplitChrPosByBin -i <input> -o output\n" \
            "Description: Convert VCF to simple code for homozyous SNP(1), heterozygous SNP(0.5), no SNP(0), missing(-1) \n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2018-08-15"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile", default=None,
                      help="Input file, CGmap or ATCGmap foramt, "
                           "use STDIN when not specified."
                           "(gzipped if end with \'gz\').", metavar="FILE")
    parser.add_option("-c", dest="codelist", default="-1,0,0.5,1",
                      help="Code list for no SNP, heterozygous SNP, homozygous SNP, seperated by \",\""
                           "Ex. a,b,c,d or 0,0,0,2. [default: %default]", metavar="STRING")
    parser.add_option("-o", dest="outfile", default=None,
                      help="Output file, use STDOUT if not specified", metavar="STRING")

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
    if (options.infile is not None) :
        try:
            if options.infile.endswith('.gz') :
                IN = gzip.open(options.infile, 'rb')
            else :
                IN = open(options.infile, 'r')
        except IOError:
            sys.stderr.write("[Error] File %s cannot be open." % fn)
            exit(-1)
    else :
        IN = sys.stdin
    #

    VCF2Matrix(IN, options.codelist)
    #
#
# ===========================================
if __name__ == "__main__":
    main()
#