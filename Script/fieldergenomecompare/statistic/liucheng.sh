#sed '1d' ABD_similar_2.csv|sed "s/\"//g"|sed "s/)//g"|sed "s/c(//g"|sed "s/,//g"|sed "s/ /\t/g" > ABD_similar_2_change.csv

#sed '1d' ABD_similar_3.csv|sed "s/\"//g"|sed "s/)//g"|sed "s/c(//g"|sed "s/,//g"|sed "s/ /\t/g" > ABD_similar_3_change.csv
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#cat `grep -v 'charac'  ABD_similar_2_change.csv | awk  '{for(i=2;i<=NF;i++) print "dev/"chr"_"$1"_"$i"_homo_level"}' chr="$chr" |tr '\n' ' '|sed s/,$//g` > ${chr}_similar_up_2 &
#cat `grep -v 'charac'  ABD_similar_3_change.csv | awk  '{for(i=2;i<=NF;i++) print "dev/"chr"_"$1"_"$i"_homo_level"}' chr="$chr" |tr '\n' ' '|sed s/,$//g` > ${chr}_similar_up_3 &
#cp `grep -v 'charac'  similar_distri/ABD_similar_1_change.csv | awk  '{for(i=2;i<=NF;i++) print "dev/"chr"_"$1"_"$i"_homo_level"}' chr="$chr" `  similar_distri/length/

#done
#sed '1d' ABD_similar_1.csv | xargs echo -n|sed "s/) /\n/g" |sed "s/,//g"|sed "s/ /\t/g"|sed "s/c(//g" >ABD_similar_1_change.csv 
#awk  '{for(i=2;i<=NF;i++) print $1"_"$i}' ABD_similar_1_change.csv|less 
#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D;do
#cat `awk '{print $1"/"CHR".homo_hmmed_snp_level"}' CHR="$CHR"  sample_path.txt | tr '\n' ' '`> ${CHR}.CN_combine_hmmed_level &
#wait_all
#done
parallel --xapply -j procfile ./bisample_compare_distance.sh ::: $(eval cat GSR_MCL/sample_path.txt) ::: ../metadata_cultivar_196_nongke.txt ::: ~/mapping/fieldergenomecompare/20200424_201_sample_compare
if false;then
for num in `seq 1 1 5`;do
num_add=$((18000001 + $num*1000000))
num_del=$((18000001 - $num*1000000))
./bisample_compare_gene_distance.sh -m ../metadata_cultivar_final.txt -v +${num}M -c ../metadata_cultivar_final.txt -h chr4D -s $num_add -b 1000000 -e $num_add -u only_homo_snp_hmm -p /data/user/shinyug/HMM_for_yzz_comp/200802/201_sample_compare &
./bisample_compare_gene_distance.sh -m ../metadata_cultivar_final.txt -v -${num}M -c ../metadata_cultivar_final.txt -h chr4D -s $num_del -b 1000000 -e $num_del -u only_homo_snp_hmm -p /data/user/shinyug/HMM_for_yzz_comp/200802/201_sample_compare &
echo $num
wait_all
done
wait
for num in `seq 1 1 5`;do
paste `awk -vn=$num '{print "Norin10_RHTB1/chr4D_"$1"_only_homo_snp_hmm-"n"M_dist"}'  ../metadata_cultivar_final.txt |tr '\n' ' '|sed s/,$//g` > Norin10_RHTB1/chr4D_rht1_-${num}_hmm_distance
paste `awk -vn=$num '{print "Norin10_RHTB1/chr4D_"$1"_only_homo_snp_hmm+"n"M_dist"}'  ../metadata_cultivar_final.txt |tr '\n' ' '|sed s/,$//g` > Norin10_RHTB1/chr4D_rht1_+${num}_hmm_distance
done
fi
