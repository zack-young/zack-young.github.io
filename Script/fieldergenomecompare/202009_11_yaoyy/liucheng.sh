#!/bin/bash
set -x
#for i in {1..7};do
#    for a in {A,B,D};do
#        paste `sed '1d' metadata_cultivar_final.txt | awk '{print $1"/chr""'"$i"'""'"$a"'"".10M.norm"}' |tr '\n' ' '|sed s/,$//g`|sed '1d' | awk '{print (NR-1)*10000000"\t"$0}' > chr$i$a.10M.join_norm
#    done
#done

#paste `awk '{print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/"$1"/combine_10M.norm"}' metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > cultivar_combine_10M_norm
#
#sed '1d' metadata_cultivar_final.txt |awk '{print $1}' |tr '\n' '\t'|sed s/,$//g |awk '{print "loc""\t"$0}' > fielder_header
#
#for i in {1..7};do     
#    for a in {A,B,D};do
#        cat fielder_header chr$i$a.1M.join_norm | sponge chr$i$a.1M.join_norm
#    done
#done
#
#for i in {1..7};do
#    for a in {A,B,D};do
#        Rscript /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/CNV_heatmap.R -i chr$i$a.1M.join_norm -M metadata_cultivar_final.txt
#    done
#done

while read line;do 
#    #a=`echo $line | awk '{print $1}'`; cd $a ; paste chr*.chr > ${a}.chr;cd ..
    a=`echo $line | awk '{print $1}'`
    sh script_10M.sh $a
#    cat $a/chr*.mask_CNV_deletion_split > $a/combine_mask_CNV_deletion_split &
#    cat $a/chr*.mask_CNV_duplication_split > $a/combine_mask_CNV_duplication_split &
#    cat $a/chr*.mask_CNV_deletion_merged > $a/combine_mask_CNV_deletion_merged &
#    cat $a/chr*.mask_CNV_duplication_merged > $a/combine_mask_CNV_duplication_merged &
#    wait_all
done < ../metadata_cultivar_final.txt
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_deletion_split"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_deletion_split
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_duplication_split"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_duplication_split
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_deletion_merged"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_deletion_merged
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_duplication_merged"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_duplication_merged

#

#
#while read line;do
#    a=`echo $line | awk '{print $1}'`
#    Rscript /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/CNV_heat_map_bySample.R -i ${a}/${a}.chr -k /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/CS_karyotype.txt -p pdf/${a}_DPbySample.pdf
#done < rest_cultivar


#while read line;do
#    for num in {1..7}; do
#        for i in {A,B,D}; do
#            for item in {deletion,duplication};do
#                a=`echo $line | awk '{print $1}'`
#                #bedtools merge -i ${a}/chr${num}${i}.mask_CNV_${item} -d 1000000 -c 4 -o distinct  > ${a}/chr${num}${i}.mask_CNV_${item}_merged &
#                ./count_CNV.py -i ${a}/chr${num}${i}.mask_CNV_${item}_merged --split on >  ${a}/chr${num}${i}.mask_CNV_${item}_split &
#            done
#        wait_all
#        done
#    done
#done < cultivar_list
paste ../bed_file/CS_4D_win1k_headed.bed `awk '{print $1"_2020-04-17162314/chr4D.1k.norm"}'  ../metadata_cultivar_final.txt` | grep -A 3000 -w '17000001' > chr4D_rht1_DP_17_20
if false;then
for num in {1..7}; do
    for i in {A,B,D}; do
        for item in {deletion,duplication}; do
            cat `awk -vnum=$num -vi=$i -vitem=$item '{print "/data/user/yangzz/mapping/fieldergenomecompare/20200416/"$1"_2020-04-17162314/chr"num""i".mask_CNV_"item}' /data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_final.txt |tr '\n' ' '|sed s/,$//g`|sort -nk1,1 - > CNV_count/chr${num}${i}.combine_mask_CNV_${item}
            ./count_CNV.py --count on --pop_num 201 -i CNV_count/chr${num}${i}.combine_mask_CNV_${item} |sort -nk1,1 >  CNV_count/chr${num}${i}.combine_mask_CNV_${item}_count 
            ./count_CNV.py --compensent on --chrom ${num}${i} -i CNV_count/chr${num}${i}.combine_mask_CNV_${item}_count -o CNV_count/chr${num}${i}.combine_mask_CNV_${item}_compensent 
            echo ${num}${i}${item}
        done
    done
done
fi



