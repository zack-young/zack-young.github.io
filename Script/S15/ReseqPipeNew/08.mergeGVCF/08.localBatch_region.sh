#!/usr/bin/env bash
set -euxo pipefail

# 2630  sh 08.localBatch_region.sh
# 2631  history
# 2632  bcftools view chr4D.raw.vcf.gz -Ob > chr4D.raw.vcf.gz
# 2633  bcftools view chr4D.raw.vcf.gz
# 2634  bcftools index chr4D.1.raw.vcf.gz
# 2635  less chr4D.1.raw.vcf.gz
# 2636  sed  's/chr4D.1/chr4D/g' chr4D.1.raw.vcf.gz >  chr4D.raw.vcf.gz
# 2637  bcftools view chr4D.raw.vcf.gz -Ob | sponge chr4D.raw.vcf.gz
# 2651  bcftools view chr4D.raw.vcf.gz --threads 4 | java -jar  /home/wangzh/bin/snpEff_latest_core/snpEff/snpEff.jar -no-downstream -no-upstream  IWGSCv1p1 - -csvStats chr4D.ann.vcf.csv | bcftools view --threads 4 -o chr4D.ann.bcf.gz -Ob
# 2652  bcftools index chr4D.ann.bcf.gz
# 2653  sh script_gene.sh -a metadata_cultivar_rht4D.txt -g TraesCS4D02G040400 use_samplelist_rht4D

SMlst=Zang1817,YM07,YM05,YM03,YM02,XM01,S97,S94,S91,S83,S78,S75,S67,S6554,S64,S6,S3331,S31,S254,S253,S252,S251,S250,S249,S248,S247,S246,S245,S244,S243,S242,S241,S240,S24,S239,S238,S237,S235,S234,S233,S232,S231,S227,S225,S224,S223,S221,S219,S218,S217,S216,S215,S214,S213,S212,S211,S210,S209,S208,S207,S206,S205,S204,S203,S188,S187,S186,S185,S184,S183,S182,S181,S180,S176,S175,S172,S170,S17,S147,S143,S142,S141,S140,S139,S137,S136,S135,S134,S133,S132,S131,S130,S129,S128,S127,S126,S125,S124,S123,S122,S12,S114,S113,S112,S111,S110,S11,S109,S108,S106,S105,S104,S103,S102,S101,S10,PH95,PH68,PH56,PH52,PH51,PH49,PH46,PH153,PH151,PH148,PH135,PH133,PH132,PH126,PH108,PH105,PH10,PH09,NZ10,NZ05,CS,C9,C8,C7,C62,C61,C60,C6,C59,C58,C57,C56,C55,C54,C53,C52,C51,C50,C5,C49,C48,C47,C46,C45,C44,C43,C41,C40,C4,C39,C38,C37,C36,C35,C34,C33,C32,C31,C30,C3,C29,C28,C27,C26,C25,C24,C23,C22,C21,C2,C19,C18,C17,C16,C15,C14,C13,C12,C11,C10,SRR10766582,SRR10766583,SRR10766588,SRR10766624,SRR10766512,SRR10766529,SRR10766564,SRR10766628,SRR10766631,ArinaLrFor,Jagger,Lancer,CDC_Landmark,Mace,Norin61,CDC_Stanley,SY_Mattis


#SMlst=HZ01,YM03
# 08
sh ./mergeGVCF_local_region.sh chr4B.1 $SMlst 30858333 30866642 #.1/2.raw.vcf.gz

# 09 & 10

#wait
# do some cleaning mannully
