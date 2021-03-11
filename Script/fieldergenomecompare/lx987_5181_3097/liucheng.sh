for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#paste `awk '{print "distance/"chr"_"$1"_torestdis"}' chr="$chr" ~/mapping/09.filterVCF/GT_AD/field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > ${chr}_119_distance
./maternal_patenal_define.py --son_define on  --i1 NZ05_S140/${chr}.homo_snp_level --i2 YM03_S140/${chr}.homo_snp_level --s1 NZ05 --s2 YM03 -o ${chr}.son_belong
for item in {NZ05,YM03,both,unknown};do

sed '1d'  ${chr}.son_belong|cut -f1,2,3,6|grep ${item} |sort -n -k 2 | bedtools merge -i -  -c 4 -o distinct > ${item}.txt

#./maternal_patenal_define.py --i1 NZ05_S140/${chr}.homo_snp_level --i2 YM03_S140/${chr}.homo_snp_level --s1 NZ05 --s2 YM03 -o ${chr}.son_belong &
done
cat {NZ05,YM03,both,unknown}.txt| sort -n -k 2  > ${chr}.son_norm
./maternal_patenal_define.py --son_norm on  --i1 ${chr}.son_norm  -o ${chr}.son_normed
./count_CNV.py --split on -i ${chr}.son_normed -o ${chr}.son_normed_splited
rm ${chr}.son_normed
rm ${chr}.son_norm
cut -f4 ${chr}.son_normed_splited | sponge ${chr}.son_normed_splited
cut -f1,2,3,4,5 ${chr}.son_belong|sort -n -k 2 | paste - ${chr}.son_normed_splited | sponge ${chr}.son_normed_splited
done
#sed '1d' ABD_similar_1.csv | xargs echo -n|sed "s/) /\n/g" |sed "s/,//g"|sed "s/ /\t/g"|sed "s/c(//g" >ABD_similar_1_change.csv 
#awk  '{for(i=2;i<=NF;i++) print $1"_"$i}' ABD_similar_1_change.csv|less 
