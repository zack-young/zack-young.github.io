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
#CHRlst=chr7B.1,chr7B.2
#,chrUn.1,chrUn.2
SMlst=YM06,YM07
#SMlst=Zang1817,YM07,YM05,YM03,YM02,XM01,S97,S94,S91,S83,S78,S75,S67,S6554,S64,S6,S3331,S31,S254,S253,S252,S251,S250,S249,S248,S247,S246,S245,S244,S243,S242,S241,S240,S24,S239,S238,S237,S235,S234,S233,S232,S231,S227,S225,S224,S223,S221,S220,S219,S218,S217,S216,S215,S214,S213,S212,S211,S210,S209,S208,S207,S206,S205,S204,S203,S188,S187,S186,S185,S184,S183,S182,S181,S180,S176,S175,S172,S170,S17,S147,S143,S142,S141,S140,S139,S137,S136,S135,S134,S133,S132,S131,S130,S129,S128,S127,S126,S125,S124,S123,S122,S12,S114,S113,S112,S111,S110,S11,S109,S108,S106,S105,S104,S103,S102,S101,S10,PH95,PH68,PH56,PH52,PH51,PH49,PH46,PH153,PH151,PH148,PH135,PH133,PH132,PH126,PH108,PH105,PH10,PH09,NZ10,NZ05,CS,C9,C8,C7,C62,C61,C60,C6,C59,C58,C57,C56,C55,C54,C53,C52,C51,C50,C5,C49,C48,C47,C46,C45,C44,C43,C41,C40,C4,C39,C38,C37,C36,C35,C34,C33,C32,C31,C30,C3,C29,C28,C27,C26,C25,C24,C23,C22,C21,C2,C19,C18,C17,C16,C15,C14,C13,C12,C11,C10


#SMlst=HZ01,YM03
# 08
if false;then
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
fi
# merge parts & annotation & index
for CHR in `for i in $(echo ${CHRlst} | tr "," " "); do echo $i|cut -d. -f1; done | sort | uniq -c | gawk '$1>=2{print $2}'`;do
    (sh ./mergeTwoParts.sh ${CHR} 
    #sh ./AnnoWithSnpEff.sh ${CHR} IWGSCv1p1 
    #bcftools index ${CHR}.ann.bcf.gz
    ) &
    wait_all
done

#wait
# do some cleaning mannully
