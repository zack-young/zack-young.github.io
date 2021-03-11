#!/usr/bin/env bash
set -euxo pipefail

# needed file: 08.mergeGVCF/mergeGVCF_local.sh
# needed file: 09.filterVCF/VCFfilter_local.sh
# needed file: 10.concatVCF/mergeTwoParts.sh
# needed file: 10.concatVCF/header_length.txt
# needed file: 10.concatVCF/AnnoWithSnpEff.sh
#lx987 S271 NZ05
#nd3097 S274 YM03
#nd5181 S140
#lx99 YM05
#jm22 NZ10
#NongDa3338 S67
#JingDong6 S6554
CHRlst=chr1A.1,chr1A.2,chr1B.1,chr1B.2,chr1D.1,chr1D.2,chr2A.1,chr2A.2,chr2B.1,chr2B.2,chr2D.1,chr2D.2,chr3A.1,chr3A.2,chr3B.1,chr3B.2,chr3D.1,chr3D.2,chr4A.1,chr4A.2,chr4B.1,chr4B.2,chr4D.1,chr4D.2,chr5A.1,chr5A.2,chr5B.1,chr5B.2,chr5D.1,chr5D.2,chr6A.1,chr6A.2,chr6B.1,chr6B.2,chr6D.1,chr6D.2,chr7A.1,chr7A.2,chr7B.1,chr7B.2,chr7D.1,chr7D.2
#CHRlst=chr5D.1,chr5D.2
#,chrUn.1,chrUn.2
SMlst=NZ33,S64,NZ34,NZ18,NZ19

#SMlst=HZ01,YM03
# 08
sh ./mergeGVCF_local.sh $CHRlst $SMlst  #.1/2.raw.vcf.gz

# 09 & 10
for CHR in `echo ${CHRlst} | tr "," " "`; do
    (sh ./VCFfilter_local.sh $CHR  #.filter.final.vcf.gz
  #
    tmp=`mktemp -d /data/user/yangzz/mapping/S15/ReseqPipeNew/08.mergeGVCF/XXXXXXX`
  #
    bcftools concat ${CHR}.indel.filter.final.vcf.gz ${CHR}.snp.filter.final.vcf.gz \
    | bcftools sort -T $tmp -O b -o ${CHR}.sort.bcf.gz ) &
    wait_all
done

wait

# merge parts & annotation & index
for CHR in `for i in $(echo ${CHRlst} | tr "," " "); do echo $i|cut -d. -f1; done | sort | uniq -c | gawk '$1>=2{print $2}'`;do
    (sh ./mergeTwoParts.sh ${CHR} 
    sh ./AnnoWithSnpEff.sh ${CHR} IWGSCv1p1 
    bcftools index ${CHR}.ann.bcf.gz) &
    wait_all
done

#wait
# do some cleaning mannully
