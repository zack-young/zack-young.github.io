while read line;do
#gawk -vchr=$CHR 'ARGIND==1{a[$1"."$2]=NR} ARGIND==2{if($1"."$2 not in a){print $0 >> chr".dist_withsm.txt"}}}' 1B_1R_region.txt $i
gawk -vchr=$line 'ARGIND==1{a[$1"."$2]=NR} ARGIND==2{if($1"."$2 in a == 0){print $0 > chr"_B_filter.txt"}}' 1B_1R_region.txt ${line}_B_mask_CNV_deletion
done < list.txt
