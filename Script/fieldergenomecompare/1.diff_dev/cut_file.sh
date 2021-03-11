num_use=$1
dev_path=$2
sample_compare_path=$3
lines=$4
CHR=$5
if [[ ! -d "${sample_compare_path}/${lines}" ]]; then
    mkdir ${sample_compare_path}/${lines}
fi

cut -f${num_use} ${dev_path}/${CHR}.1M.combinediff > ${sample_compare_path}/${lines}/${CHR}.1M.density
paste -d '\t' ${CHR}.1M.bed.tmp ${sample_compare_path}/${lines}/${CHR}.1M.density | sponge ${sample_compare_path}/${lines}/${CHR}.1M.density
