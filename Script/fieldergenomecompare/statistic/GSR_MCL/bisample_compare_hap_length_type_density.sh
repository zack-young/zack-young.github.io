#!/bin/bash
#set -euxo pipefail
CHR=$1
n=$2
START=$(( ($n-1)*1000000+1 ))
sample=$3
if false;then
line=$3
SAMPLE1=`echo $line|awk -F ',' '{print $1}'`
SAMPLE2=`echo $line|awk -F ',' '{print $2}'`
START=$(( ($n-1)*1000000+1 ))
DEV="~/mapping/fieldergenomecompare/20200424_201_sample_compare/${SAMPLE1}_${SAMPLE2}"
  if [ `wc -l ${DEV}/${CHR}.hmm.v4.low_both_merge.txt` !=0 ];then
    length=`awk '-vSTART=${START} {a=0;if(START>=$2&&START<=$3){a=$3-$2}}END{print a}' ~/mapping/fieldergenomecompare/20200424_201_sample_compare/${SAMPLE1}_${SAMPLE2}/${CHR}.hmm.v4.low_both_merge.txt`
  else
    length=0
  fi
if [ ${length} != 0 ];then
echo -e "${CHR}\t${START}\t${length}" >> ${CHR}.hap_length
fi
fi

if false;then
line_CNLC=`awk -vSTART=${START} '{if ($2==START){print $0}}' ${CHR}_ALL_hap_count.txt |grep -vE 'exclusive|deletion|duplication'| grep -vE 'ASC|AFC|AMC|AUC|EUC'|grep 'CNL'| grep 'CNC'| wc -l `
line_CNCF=`awk -vSTART=${START} '{if ($2==START){print $0}}' ${CHR}_ALL_hap_count.txt |grep -vE 'exclusive|deletion|duplication'|grep -E 'ASC|AFC|AMC|AUC|EUC'| grep 'CNC' | grep -v 'CNL'| wc -l`
line_CNCA=`awk -vSTART=${START} '{if ($2==START){print $0}}' ${CHR}_ALL_hap_count.txt |grep -vE 'exclusive|deletion|duplication'|grep -E 'ASC|AFC|AMC|AUC|EUC' | grep 'CNC' | grep 'CNL'| wc -l`

line1_CNL=`awk -vSTART=${START} '{if ($2==START){print $0}}' ${CHR}_CNL_hap_count.txt | grep -vE 'exclusive|deletion|duplication' | wc -l`
line2_CNL=`awk -vSTART=${START} '{if ($2==START){print $0}}' ${CHR}_CNL_hap_count.txt | grep 'exclusive' | cut -f5 | tr ',' '\n' | wc -l`

line1_CNC=`awk -vSTART=${START} '{if ($2==START){print $0}}' ${CHR}_CNC_hap_count.txt | grep -vE 'exclusive|deletion|duplication' | wc -l`
line2_CNC=`awk -vSTART=${START} '{if ($2==START){print $0}}' ${CHR}_CNC_hap_count.txt | grep 'exclusive' | cut -f5 | tr ',' '\n' | wc -l`


echo ${line1_CNC}	${line2_CNC}	${line1_CNL}	${line2_CNL}	${line_CNLC}	${line_CNCF} ${line_CNCA}| awk -vCHR=$CHR -vSTART=$START '{print CHR"\t"START"\t"$1+$2"\t"$3+$4"\t"$5"\t"$6"\t"$7}' >> CNC_CNL_ALL_hap_stat_density.txt
fi
#===========
if true;then
line_CNL=`cat chr*_ALL_hap_count.txt |grep -vE 'exclusive|deletion|duplication'|grep -vE 'ASC|AFC|AMC|AUC|EUC|ASL|AFL|AML|AUL|EUL|CNW'|grep 'CNL'|grep $sample |wc -l`
line_F=`cat chr*_ALL_hap_count.txt |grep -vE 'exclusive|deletion|duplication'|grep -vE 'CNL|CNW'|grep -E 'ASC|AFC|AMC|AUC|EUC|ASL|AFL|AML|AUL|EUL'|grep $sample |wc -l`
line_CNC1=`cat chr*_ALL_hap_count.txt |grep -vE 'exclusive|deletion|duplication'|grep -vE 'ASC|AFC|AMC|AUC|EUC|ASL|AFL|AML|AUL|EUL|CNW|CNL'|grep $sample |wc -l`
line_CNC2=`cat chr*_ALL_hap_count.txt |grep 'exclusive'|grep $sample|wc -l`
#line_BOTH=`cat chr*_ALL_hap_count.txt |grep -vE 'exclusive|deletion|duplication'| grep -E 'ASC|AFC|AMC|AUC|EUC'|grep 'CNL'|grep $sample |wc -l`
#line_FL=`cat chr*_ALL_hap_count.txt |grep -vE 'exclusive|deletion|duplication'|grep -E 'ASL|AFL|AML|AUL|EUL'|grep -vE 'CNL|ASC|AFC|AMC|AUC|EUC|CNW'|grep $sample |wc -l`
line_FC_CNL=`cat chr*_ALL_hap_count.txt |grep -vE 'exclusive|deletion|duplication'|grep -vE 'ASL|AFL|AML|AUL|EUL|CNW'| grep -E 'ASC|AFC|AMC|AUC|EUC'|grep CNL |grep $sample |wc -l`
#line_CNW=`cat chr*_ALL_hap_count.txt |grep -vE 'exclusive|deletion|duplication'| grep -vE 'ASC|AFC|AMC|AUC|EUC|ASL|AFL|AML|AUL|EUL|CNL'|grep 'CNW'|grep $sample |wc -l`
line_CNC=`echo $line_CNC1	$line_CNC2 | awk '{print $1+$2}'`
line_oth=`cat chr*_ALL_hap_count.txt |grep -vE 'deletion|duplication'|grep $sample |wc -l`

echo -e "${sample}\t${line_CNL}\t${line_F}\t${line_CNC}\t${line_FC_CNL}\t${line_oth}" >> CNCsample_2_CNL_F.txt
fi

#============
if false;then
line1_CNL=`awk -vSTART=${START} '{if ($2==START){print $0}}' ${CHR}_ALL_hap_count.txt | grep -vE 'exclusive|deletion|duplication' | wc -l`
line2_CNL=`awk -vSTART=${START} '{if ($2==START){print $0}}' ${CHR}_ALL_hap_count.txt | grep 'exclusive' | cut -f5 | tr ',' '\n' | wc -l`
echo ${line1_CNL}	${line2_CNL} | awk -vCHR=$CHR -vSTART=$START '{print CHR"\t"START"\t"$1+$2}' >> ALL_hap_stat_density.txt
fi

