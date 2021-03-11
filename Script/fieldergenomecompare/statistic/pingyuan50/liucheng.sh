sed  's/[^\t]\{1,5\}duplication_CNV\>/10/g' 5D_heatmap|sed  's/[^\t]\{1,5\}deletion_CNV\>/-10/g'|less -S

$ sed  's/S10duplication_CNV\>/7/g' 5D_heatmap|sed  's/S10deletion_CNV\>/8/g' >  5D_heatmap_change 

$ sed  's/duplication_both_CNV/4/g' 5D_heatmap_change|sed  's/deletion_both_CNV\>/5/g'|sponge 5D_heatmap_change

$ sed  's/[^\t]\{1,5\}duplication_CNV\>/10/g' 5D_heatmap_change|sed  's/[^\t]\{1,5\}deletion_CNV\>/9/g'|sponge 5D_heatmap_change

$ sed  's/low/0/g' 5D_heatmap_change|sed  's/mid/1/g'| sed  's/high/2/g'|sponge 5D_heatmap_change

while read line; do
sort -n -k 2 ~/mapping/fieldergenomecompare/20200424_201_sample_compare/${line}/chr5D.homo_snp_level|cut -f 4 > chr5D.homo_snp_level_${line}
done < sample_list

paste `awk '{print "chr5D.homo_snp_level_"$1}'  sample_list |tr '\n' ' '|sed s/,$//g` > 5D_heatmap
cat sample_list| tr '\n' '\t'  |sed s/\\t$/\\n/g |cat - 5D_heatmap_change|sponge 5D_heatmap_change

sed  's/S10_//g' 5D_heatmap|sed  's/_S10//g' |sponge 5D_heatmap
