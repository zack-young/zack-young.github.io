num=1
for dev in CP06_CP08 CP06_CP11 CP06_884187 CP07_CP09 CP11_CP09 CP09_884187; do
if [ "$num" == 4 ];then
num=1
fi
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
sort -k 2,2 -n ${dev}/${CHR}.homo_undefined_snp_level | sponge ${dev}/${CHR}.homo_undefined_snp_level
sed s/'low'/"low${num}"/g ${dev}/${CHR}.homo_undefined_snp_level > ${dev}/${CHR}.snp_level_sp
done
num=`expr $num + 1`
done

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
paste CP06_CP08/${CHR}.homo_undefined_snp_level  CP06_CP11/${CHR}.homo_undefined_snp_level CP06_884187/${CHR}.homo_undefined_snp_level | cut -f1,2,3,4,8|awk '{if ($0!~"low"&& $0!~"both" ){print $1"\t"$2"\t"$3"\t""undefined"}}' > CP06_CP08_CP11_884187/${CHR}.undefined_level
paste CP07_CP09/${CHR}.homo_undefined_snp_level  CP11_CP09/${CHR}.homo_undefined_snp_level CP09_884187/${CHR}.homo_undefined_snp_level | cut -f1,2,3,4,8|awk '{if ($0!~"low"&& $0!~"both" ){print $1"\t"$2"\t"$3"\t""undefined"}}' > CP07_CP09_CP11_884187/${CHR}.undefined_level
done

../statistic/NongDa5181/ptf_maker.py -p '/data/user/yangzz/mapping/fieldergenomecompare/202009_11_yaoyy/CP06_CP08_CP11_884187,/data/user/yangzz/mapping/fieldergenomecompare/202009_11_yaoyy/CP06_CP08,/data/user/yangzz/mapping/fieldergenomecompare/202009_11_yaoyy/CP06_CP11,/data/user/yangzz/mapping/fieldergenomecompare/202009_11_yaoyy/CP06_884187' -n 4 -s ".undefined_level,.snp_level_sp,.snp_level_sp,.snp_level_sp" > plotfile_CP06_CP08_CP11_884
