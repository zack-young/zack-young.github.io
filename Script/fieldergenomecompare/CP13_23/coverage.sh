i=$1
while read line;do
echo $line|awk  -F " " '{print $1"\t"$2"\t"$3"\t"$4}' > test${i}.txt
chr=`echo $line|awk -F " " '{print $1}'`
a=`bedtools coverage -a test${i}.txt -b ~/mapping/Reseq_data/CP${i}/05.mergeAsplit/${chr}.*.bam -counts -sorted`
echo $a 'CP'$i >> gene_use.txt
done < gene2.txt
