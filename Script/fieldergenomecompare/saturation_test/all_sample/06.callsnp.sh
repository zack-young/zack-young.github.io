#!/usr/bin/env bash
# Guo, Weilong; guoweilong@126.com; 2017-10-24
#echo $CHR
chrlst=chr1A.1,chr1A.2,chr1B.1,chr1B.2,chr1D.1,chr1D.2,chr2A.1,chr2A.2,chr2B.1,chr2B.2,chr2D.1,chr2D.2,chr3A.1,chr3A.2,chr3B.1,chr3B.2,chr3D.1,chr3D.2,chr4A.1,chr4A.2,chr4B.1,chr4B.2,chr4D.1,chr4D.2,chr5A.1,chr5A.2,chr5B.1,chr5B.2,chr5D.1,chr5D.2,chr6A.1,chr6A.2,chr6B.1,chr6B.2,chr6D.1,chr6D.2,chr7A.1,chr7A.2,chr7B.1,chr7B.2,chr7D.1,chr7D.2
DIR='/data2/user2/xiexm/pro_ject/resequencing/compare/s30/05.mergeAsplit'
DIR1='all_sample'
CHR=`echo ${chrlst} | tr "," " "`
parallel -j procfile --results ${DIR1} sh csp.sh ::: $(eval echo ${CHR}) ::: ${DIR} 
