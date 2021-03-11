#!/usr/bin/env python

def Cluster ( infile="inflie.txt", bin_size = 1000000, cutoff = 2000, DP_file="DPmatrix.txt", GQ_file='GQmatrix.txt', outfile="outfile.txt"):

    x=np.loadtxt(infile, "i2")  # int8, int16, int32, int64 四种数据类型可以使用字符串 'i1', 'i2','i4','i8' 代替
    DP_matrix=np.loadtxt(DP_file, "i2")
    GQ_matrix=np.loadtxt(GQ_file, "i2")
    DP_final=np.where((DP_matrix<3)|(DP_matrix>99),30,1)
    GQ_final=np.where(GQ_matrix<8,30,1)
    final=x*GQ_final
    final=final*DP_final
    s_result=abs(final[:,None,:]-final[:,:,None])
    s_result[s_result>5]=0
    z=s_result.sum(axis=0)
    array=z[np.triu_indices(z.shape[0],k=1)]/4
    ## x[:,:,None] equal to numpy.newaxis and by this way we can make sure subtract between items
    np.savetxt(outfile, array, newline='\t', fmt="%.8f")
    #geno_mat=np.loadtxt(infile+".rawgeno", "i2")
    #missing_mat=np.zeros(geno_mat.shape, "i2")  #return an array filled by 0 in specified shape
    #missing_mat[geno_mat<3]=1
    #missing_rate=sum(missing_mat)/geno_mat.shape[0]
    #geno_mat[geno_mat>2]=0
    #z=sum(geno_mat)/missing_rate
    #z[missing_rate==0]=0

    #geno_mat_pair=np.loadtxt(DIST_PATH+infile+".rawdist.dist")
    #geno_mat_pair=np.hstack((np.transpose([z]), geno_mat_pair))
    #geno_mat_pair=np.concatenate(([np.append(0,z)],geno_mat_pair))
    #geno_mat_pair=np.nan_to_num(geno_mat_pair)/2
    #np.savetxt(infile+".dist_to_CS.txt", geno_mat_pair, "%i")



from optparse import OptionParser
import numpy as np
from sklearn.cluster import AgglomerativeClustering

# ===========================================
def main():
    parser = OptionParser()
    parser.add_option("-i", dest="infile")
    parser.add_option("-b", dest="bin_size", default=1000000)
    parser.add_option("-c", dest="cutoff", default=2000)
    parser.add_option("-o", dest="outfile")
    parser.add_option("-d", dest="DP_file")
    parser.add_option("-g", dest="GQ_file")
    (options, args) = parser.parse_args()
    #
    Cluster (options.infile, int(options.bin_size), int(options.cutoff), options.DP_file, options.GQ_file, options.outfile)

# ===========================================
if __name__ == "__main__":
    main()
