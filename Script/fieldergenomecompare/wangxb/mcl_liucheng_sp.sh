BED=$1
for CHR in `cut  -f1 ${BED}|sort | uniq`; do
#CHR="chr1A"
n=`grep $CHR $BED | wc -l`
start=`grep $CHR $BED|cut -f2 | tr '\n' ' '`
#parallel --xapply -j procfile ./bisample_compare_bin_edge_list.sh sample_path.txt ::: $(eval echo {1..$n}) ::: ${CHR} ::: $start ::: homo_undefined_snp_level ::: /data/user/yangzz/mapping/fieldergenomecompare/wangxb/dev ::: statistic ::: meta_data_all.txt
for i in $(eval echo {1..$n});do
/home/wangzh/bin/mcl statistic/${CHR}_homo_undefined_snp_level_${i}  --abc -o statistic/${CHR}_mcl_${i}
done
done
#./bisample_compare_bin_edge_list.sh  sample_path.txt  ppd chr2D 33500001  homo_undefined_snp_level ~/mapping/fieldergenomecompare/ppd1  ppd_mcl ../metadata_cultivar_196_nongke_10+.txt
