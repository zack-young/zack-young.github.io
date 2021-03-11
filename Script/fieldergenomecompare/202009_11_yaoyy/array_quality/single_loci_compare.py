#!/data/user/yangzz/worktools/anaconda3/bin/python
# #!/usr/bin/env python
from scipy import stats
def Cluster ( infile="inflie.txt", bin_size = 1000000, head_file = "file", DP_file="DPmatrix.txt", GQ_file='GQmatrix.txt',gene_sample_type='sample_list',num=1 ,CHR="chr1A",outfile="outfile.txt"):
    #from itertools import combinations
    #for i in list(combinations(typ_sample, 2)):
    sample_lis=[]
    with open(head_file) as f:
        sample=f.readlines()
    for i in sample:
        tokens=i.strip().split("\t")
        sample_lis.append(tokens[2])
###------------------------------------------###
    typ_sample=gene_sample_type.split("|")
    for item_all in typ_sample:
        sample_index=[]
        tokens=item_all.split(",")
        for item in tokens:
            sample_index.append(sample_lis.index(item[4:]))

    
    #CNV_PATH="/data/user/yangzz/mapping/fieldergenomecompare/5.cross_sample/CNV_masker/"
        x=np.loadtxt(infile, "i2")  # int8, int16, int32, int64 四种数据类型可以使用字符串 'i1', 'i2','i4','i8' 代替
        #cnv_mat=np.loadtxt(CNV_file)
        #cnv_mat_tmp=cnv_mat[num-1,]
        DP_matrix=np.loadtxt(DP_file, "i2")
        GQ_matrix=np.loadtxt(GQ_file, "i2")
        DP_final=np.where((DP_matrix<3)|(DP_matrix>99),30,1)
        GQ_final=np.where(GQ_matrix<8,30,1)
        missing_mat=np.zeros(x.shape, "i2")
        missing_mat[(DP_final<30)&(GQ_final<30)&(x<30)]=1
        #print(len(missing_mat))
        final=x*DP_final
        final=final*GQ_final
        use_array=final[sample_index]
        final_filter=np.delete(final,sample_index, axis=0)
        fail=0
        sim=0
        diff=0
        s_result=abs(use_array[:,None]-use_array[None,:])
        s_result[s_result>5]=4
        array=s_result[np.triu_indices(s_result.shape[0],k=0)]/4
        if len(array)==0:
            sim_ratio=0
        else:
            sim_ratio=sum(array)/len(array) 
        mode_num=stats.mode(use_array)[0][0]
        if 1-sim_ratio>=0.9 and len(use_array[use_array<11])>=len(sample_index)*0.9 :
            for i in final_filter:
                if i >= 30:
                     fail += 1
                elif i == mode_num:
                     sim += 1
                else:
                     diff += 1
            ratio_diff=diff/(fail+sim+diff)
            ratio_fail=fail/(fail+sim+diff)
            print('\t'.join([CHR,str(num),",".join(tokens),"%.3f" % sim_ratio,"%.3f" % ratio_diff,"%.3f" % ratio_fail]))


from optparse import OptionParser
import numpy as np
from sklearn.cluster import AgglomerativeClustering

# ===========================================
def main():
    parser = OptionParser()
    parser.add_option("-i", dest="infile")
    parser.add_option("-b", dest="bin_size", default=1000000)
    parser.add_option("-u", dest="head_file", default="file")
    parser.add_option("-e", dest="gene_sample_type")
    parser.add_option("-c", dest="CHR")
    parser.add_option("-n", dest="num")
    parser.add_option("-o", dest="outfile")
    parser.add_option("-d", dest="DP_file")
    parser.add_option("-g", dest="GQ_file")
    (options, args) = parser.parse_args()
    #
    Cluster (options.infile, int(options.bin_size), options.head_file, options.DP_file, options.GQ_file, options.gene_sample_type, int(options.num),options.CHR, options.outfile)

# ===========================================
if __name__ == "__main__":
    main()
