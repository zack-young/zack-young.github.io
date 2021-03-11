total_num=`bcftools query -l chr1B.S.rawgeno| wc -l `
for num in `seq 1 $total_num`;do
(sample_num=`cut -f${num} chr1B.S.matrix | grep '30' |wc -l`
sample=`cut -f${num} header.txt`
final_num=`awk 'BEGIN{printf "%.2f\n",('$sample_num'/192852)}'`
echo ${sample}	${final_num} >> chr1B.S.miss_count) &
sleep 0.001s
done
