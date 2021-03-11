#!/bin/bash
#set -x
if false;then
for i in {1..7};do
    for a in {A,B,D};do
        paste ` awk '{print "/data/user/yangzz/mapping/fieldergenomecompare/20200416_201_CNV/"$1"_2020-04-17162314/chr""'"$i"'""'"$a"'"".1M.norm_un1k"}' ../metadata_cultivar_196_nongke_10+.txt |tr '\n' ' '|sed s/,$//g`|sed '1d' | awk '{print (NR-1)*1000000"\t"$0}' > chr$i$a.1M.join_norm_un1k
    done
done

#paste `awk '{print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/"$1"/combine_10M.norm"}' metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > cultivar_combine_10M_norm
#
awk '{print $1}' ../metadata_cultivar_196_nongke_10+.txt |tr '\n' '\t'|sed s/"\t"$//g |awk '{print "loc""\t"$0}' > 196_header
#
for i in {1..7};do     
    for a in {A,B,D};do
        cat 196_header chr$i$a.1M.join_norm_un1k | sponge chr$i$a.1M.join_norm_un1k
    done
done
cat chr*.1M.join_norm_un1k  | grep -v 'loc'|awk -F "\t"  '{for (i=2;i<=NF-1;i++)printf("%s\t", $i);print $NF}'> combine.1M.join_norm_un1k
awk -F "\t"  '{for (i=2;i<=NF-1;i++)printf("%s\t", $i);print $NF}' 196_header |cat - combine.1M.join_norm | sponge combine.1M.join_norm
fi

if false;then
for item in {duplication,deletion};do
for i in `awk '{print $1}' ~/mapping/fieldergenomecompare/metadata_cultivar_196_nongke_10+.txt`;do
cat /data/user/yangzz/mapping/fieldergenomecompare/20200416_201_CNV/${i}_2020-04-17162314/chr*A.mask_CNV_${item} > subgenome_CNV_count/${i}_A_mask_CNV_${item}
cat /data/user/yangzz/mapping/fieldergenomecompare/20200416_201_CNV/${i}_2020-04-17162314/chr*B.mask_CNV_${item} > subgenome_CNV_count/${i}_B_mask_CNV_${item}
cat /data/user/yangzz/mapping/fieldergenomecompare/20200416_201_CNV/${i}_2020-04-17162314/chr*D.mask_CNV_${item} > subgenome_CNV_count/${i}_D_mask_CNV_${item}
done
wc -l subgenome_CNV_count/*A_mask_CNV_${item} |grep "mask_CNV_${item}" > A_${item}_number
wc -l subgenome_CNV_count/*B_mask_CNV_${item} |grep "mask_CNV_${item}" > B_${item}_number
wc -l subgenome_CNV_count/*D_mask_CNV_${item} |grep "mask_CNV_${item}" > D_${item}_number
done
fi
#wc -l *A_mask_CNV_deletion |grep 'mask_CNV_deletion' > A_deletion_number
#
#for i in {1..7};do
#    for a in {A,B,D};do
#        Rscript /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/CNV_heatmap.R -i chr$i$a.1M.join_norm -M metadata_cultivar_final.txt
#    done
#done

#while read line;do 
#    #a=`echo $line | awk '{print $1}'`; cd $a ; paste chr*.chr > ${a}.chr;cd ..
#    a=`echo $line | awk '{print $1}'`
#    cat $a/chr*.mask_CNV_deletion_split > $a/combine_mask_CNV_deletion_split &
#    cat $a/chr*.mask_CNV_duplication_split > $a/combine_mask_CNV_duplication_split &
#    cat $a/chr*.mask_CNV_deletion_merged > $a/combine_mask_CNV_deletion_merged &
#    cat $a/chr*.mask_CNV_duplication_merged > $a/combine_mask_CNV_duplication_merged &
#    wait_all
#done < cultivar_list
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

if true;then
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
        for item in {deletion,duplication}; do
            for typ in {CN,CNC,CNL,ALL};do
            #cat `awk -vnum=$num -vi=$i -vitem=$item '{print "/data/user/yangzz/mapping/fieldergenomecompare/20200416_201_CNV/"$1"_2020-04-17162314/chr"num""i".mask_CNV_"item}' /data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_${typ}.txt |tr '\n' ' '|sed s/,$//g`|sort -nk1,1 - > CNV_${typ}_count/chr${num}${i}.combine_mask_CNV_${item}
            #line_num=`wc -l /data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_${typ}.txt`
            cat `awk -vCHR=$CHR  -vitem=$item '{print "/data/user/yangzz/mapping/fieldergenomecompare/20200416_201_CNV/"$1"_2020-04-17162314/"CHR".mask_CNV_"item}' /data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_${typ}.txt |tr '\n' ' '|sed s/,$//g`|sort -nk2,2 - > CNV_${typ}_count/${CHR}.combine_mask_CNV_${item}
            line_num=`wc -l /data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_196_nongke_10+.txt`
            cut -f2,3 CNV_${typ}_count/${CHR}.combine_mask_CNV_${item} | uniq -c | awk '{print $2"\t"$3"\t"$1}'> CNV_${typ}_count/${CHR}.combine_mask_CNV_${item}_count
            count=`cat CNV_${typ}_count/${CHR}.combine_mask_CNV_${item}|wc -l `
            echo $count
            if [ $count != 0 ];then
              while read line ;do
                a=`echo $line|awk '{print $1}'`
                b=`awk -va=$a '{if($2==a) print $5}' CNV_${typ}_count/${CHR}.combine_mask_CNV_${item}| tr '\n' ':'|sed s/:$//g`
                echo $line $b >> CNV_${typ}_count/${CHR}.combine_mask_CNV_${item}_count_sample
              done<CNV_${typ}_count/${CHR}.combine_mask_CNV_${item}_count
            else
              touch CNV_${typ}_count/${CHR}.combine_mask_CNV_${item}_count_sample
            fi
            #./count_CNV.py --ratio on --pop_num ${line_num} -i CNV_${typ}_count/chr${num}${i}.combine_mask_CNV_${item} |sort -nk1,1 >  CNV_${typ}_count/chr${num}${i}.combine_mask_CNV_${item}_ratio
            #./count_CNV.py --count on --pop_num ${line_num} -i CNV_${typ}_count/chr${num}${i}.combine_mask_CNV_${item} |sort -nk1,1 >  CNV_${typ}_count/chr${num}${i}.combine_mask_CNV_${item}_count 
            #./count_CNV.py --compensent on --chrom ${num}${i} -i CNV_${typ}_count/chr${num}${i}.combine_mask_CNV_${item}_ratio -o CNV_${typ}_count/chr${num}${i}.combine_mask_CNV_${item}_ratio_compensent 
            #./count_CNV.py --compensent on --chrom ${num}${i} -i CNV_${typ}_count/chr${num}${i}.combine_mask_CNV_${item}_count -o CNV_${typ}_count/chr${num}${i}.combine_mask_CNV_${item}_count_compensent
            echo ${CHR}${item}
            done
        done
done
fi

#for num in `awk  -vitem=$item '{print $1}' /data/user/yangzz/mapping/fieldergenomecompare/metadata_cultivar_CNC.txt |tr '\n' ' '|sed s/,$//g`;do
#awk -vitem=$num '{print $0"\t"item}' /data/user/yangzz/mapping/fieldergenomecompare/20200416_201_CNV/${num}_2020-04-17162314/chr5A.1M.norm | sed -n '589p'
#done
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
cat CNV_CNL_count/${CHR}.combine_mask_ALL_CNV_count_compensent | awk '{if($5 > 51) print chr"\t"$0}' chr="$CHR" >> noise_CNV.txt
done
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
cat CNV_CNC_count/${CHR}.combine_mask_ALL_CNV_count_compensent | awk '{if($5 > 53) print chr"\t"$0}' chr="$CHR" >> noise_CNV.txt
done
sort -k1,1 -k2,2n noise_CNV.txt | cut -f1,2,3 |uniq |sponge noise_CNV.txt
