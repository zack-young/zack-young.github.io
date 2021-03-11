#!/usr/bin/env python

from optparse import OptionParser
import numpy as np
import glob

def concat_dist(CHR="chr1A", BINSIZE=5, THR=1000, OUTFILE="test"):
    #
    BED_PATH="/data2/rawdata2/tetraintro/bed_file2/"
    DIST_PATH="/data2/rawdata2/variant_density/6.cross_sample_5M/raw/"
    init_mat=np.zeros((366, 366))
    count=0
    
    for CHR in ["chr1A","chr1B","chr2A","chr2B","chr3A","chr3B","chr4A","chr4B","chr5A","chr5B","chr6A","chr6B","chr7A","chr7B"]:
        for FILE in glob.glob(CHR + "*.txt"):
            tmp=np.loadtxt(DIST_PATH+FILE)
            init_mat=init_mat+tmp
            count+=1
    #
    init_mat=init_mat/count
    np.savetxt("SM_level.dist.txt", init_mat, "%i")


def main():
    ''' convert (concat) dist between bins to dist between samples.'''
    parser = OptionParser()
    parser.add_option("-c", dest="CHR")
    parser.add_option("-b", dest="BINSIZE", default=5)
    parser.add_option("-t", dest="THR", default=1000)
    parser.add_option("-o", dest="OUTFILE")
    (options, args) = parser.parse_args()
    #
    concat_dist(options.CHR, options.BINSIZE, options.THR, options.OUTFILE)
    
# ==========================================
if __name__ == "__main__":
    main()
