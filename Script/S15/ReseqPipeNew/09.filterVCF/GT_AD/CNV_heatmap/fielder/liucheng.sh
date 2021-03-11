for i in {1..7};do     
for a in {A,B,D};do
cat fielder_header chr$i$a.1M.join_norm | sponge chr$i$a.1M.join_norm
done
done

for i in {1..7};do
for a in {A,B,D};do
Rscript ../bin/CNV_heatmap.R -i chr$i$a.1M.join_norm -M ../../field_cultivar/metadata_cultivar.txt
done
done

for i in {1..7};do\n    for a in {A,B,D};do\npaste `awk '{print "/data2/rawdata2/readDepth/"$1"/chr""'"$i"'""'"$a"'"".1M.norm"}' metadata_cultivar.txt |tr '\n' ' '|sed s/,$//g`| sed '1d' | awk '{print (NR-1)*1000000"\t"$0}' > chr$i$a.1M.join_norm\ndone\ndone

while read line;do a=`echo $line | awk '{print $1}'`; cd $a ; paste chr*.chr > ${a}.chr;cd ..; done < metadata_cultivar.txt

while read line;do\na=`echo $line | awk '{print $1}'`\nfor chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do\nmkdir ${a} ; sed '1 s/^.*$/'$chr'/' /data2/rawdata2/readDepth/${a}/${chr}.1M.norm > ${a}/${chr}.1M.norm.chr\ndone\ndone < metadata_cultivar.txt
