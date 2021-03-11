num_use=$1
dev_path=$2
sample_compare_path=$3
lines=$4
CHR=$5
BED=$6
if [[ ! -d "${sample_compare_path}/${lines}" ]]; then
    mkdir ${sample_compare_path}/${lines}
fi

cut -f${num_use} ${dev_path}/${CHR}.1M.combinediff > ${sample_compare_path}/${lines}/${CHR}.1M.density
grep ${CHR} ${BED} > ${num_use}.bed
sed -i "1i\CHR\tstart\tend" ${num_use}.bed
cut -f1,2,3 ${num_use}.bed|paste -d '\t' - ${sample_compare_path}/${lines}/${CHR}.1M.density | sponge ${sample_compare_path}/${lines}/${CHR}.1M.density
rm -f ${num_use}.bed
