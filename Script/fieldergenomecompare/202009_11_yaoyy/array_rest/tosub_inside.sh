CHR=$1
START=$2
END=$3
N=$4
meta_file=$5
sample_lis=$6
ord=$7


tmpfile_gt=`mktemp -p /data/user/yangzz/mapping/fieldergenomecompare/statistic/heatmap/`
tmpfile_GQ=`mktemp -p /data/user/yangzz/mapping/fieldergenomecompare/statistic/heatmap/`
tmpfile_DP=`mktemp -p /data/user/yangzz/mapping/fieldergenomecompare/statistic/heatmap/`
awk -vN=${N} 'NR==1 || (NR==N+1){print $0}' rawdata/${CHR}.${START}-${END}.matrix > ${tmpfile_gt}
awk -vN=${N} 'NR==1 || (NR==N+1){print $0}' rawdata/${CHR}.${START}-${END}.DPmatrix > ${tmpfile_DP}
awk -vN=${N} 'NR==1 || (NR==N+1){print $0}' rawdata/${CHR}.${START}-${END}.GQmatrix >${tmpfile_GQ}
pos=`tail -n +${N} rawdata/${CHR}.${START}-${END}.POS|head -1`
#echo $pos
python single_loci_compare.py -i ${tmpfile_gt} -d ${tmpfile_DP} -g ${tmpfile_GQ} -u ${meta_file} -e ${sample_lis} -n ${pos} -c ${CHR} >> loci_raw${ord}.txt

rm -f ${tmpfile_gt}
rm -f ${tmpfile_GQ}
rm -f ${tmpfile_DP}
