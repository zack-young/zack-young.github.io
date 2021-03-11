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
#CHRlst=chr1A.1,chr1A.2,chr1B.1,chr1B.2,chr1D.1,chr1D.2,chr2A.1,chr2A.2,chr2B.1,chr2B.2,chr2D.1,chr2D.2,chr3A.1,chr3A.2,chr3B.1,chr3B.2,chr3D.1,chr3D.2,chr4A.1,chr4A.2,chr4B.1,chr4B.2,chr4D.1,chr4D.2,chr5A.1,chr5A.2,chr5B.1,chr5B.2,chr5D.1,chr5D.2,chr6A.1,chr6A.2,chr6B.1,chr6B.2,chr6D.1,chr6D.2,chr7A.1,chr7A.2,chr7B.1,chr7B.2,chr7D.1,chr7D.2,chrUn.1,chrUn.2
CHRlst=chr1A.1,chr1A.2,chr1B.1,chr1B.2,chr1D.1,chr1D.2,chr2A.1,chr2A.2,chr2B.1,chr2B.2,chr2D.1,chr2D.2,chr3B.1,chr3B.2,chr3D.1,chr3D.2,chr4A.1,chr4A.2,chr4B.1,chr4B.2,chr4D.1,chr4D.2,chr5A.1,chr5A.2,chr5B.1,chr5B.2,chr5D.1,chr5D.2,chr6A.1,chr6A.2,chr6B.1,chr6B.2,chr7A.1,chr7A.2,chr7B.1,chr7B.2,chr7D.1,chr7D.2,chrUn.1,chrUn.2
#SMlst=S12,S234,S235,C1,C2,C3,S233,C24,C25,C60,S128,S129,S134,S204,C5,C6,C7,C8,C9,C10,S132,S123,S131,S203,C31,C33,C32,C34,C35,S206,S207,S208,C20,C21,S122,C26,C27,C28,PH09,C11,C12,C13,C14,C15,C16,C17,C18,C19,C48,C49,C58,S231,S232,C4,S212,S130,S236,S237,S238,NZ01,NZ05,NZ10,NZ11,YM03,CS,C41,C42,C43,C44,C45,C46,C47,S13,S64,S67,S135,S136,S137,S139,S140,S142,S209,S210,S211,S213,S214,S215,S216,S217,S218,S219,S220,S221,S223,S224,S225,S227,S228,S229,S230,S3331,S6554,C29,C30,S133,S126,S127,C50,C51,C52,C53,C54,C55,C56,C57,C59,C61,C62,S239,C22,C23,S125,S205,S124
SMlst=PH16_to_Spelt,PH17_to_Spelt
# 08
sh ./mergeGVCF_local.sh $CHRlst $SMlst

# 09 & 10
for CHR in `echo ${CHRlst} | tr "," " "`; do
    (sh ./VCFfilter_local.sh $CHR
  #
    tmp=`mktemp -d /data/user/yangzz/mapping/S15/ReseqPipeNew/08.mergeGVCF/XXXXXXX`
  #
    bcftools concat ${CHR}.indel.filter.final.vcf.gz ${CHR}.snp.filter.final.vcf.gz \
    | bcftools sort -T $tmp -O b -o ${CHR}.sort.bcf.gz ) &
    wait_all
done

wait
#
# merge parts & annotation & index
for CHR in `for i in $(echo ${CHRlst} | tr "," " "); do echo $i|cut -d. -f1; done | sort | uniq -c | gawk '$1>=2{print $2}'`;do
    sh ./mergeTwoParts.sh ${CHR} 
    sh ./AnnoWithSnpEff.sh ${CHR} IWGSCv1p1 
    bcftools index ${CHR}.ann.bcf.gz
    wait_all
done

wait
# do some cleaning mannully
