sample1=${1}
sample2=${2}
son=${3}
if false;then
for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#for chr in chr4B;do
#paste `awk '{print "distance/"chr"_"$1"_torestdis"}' chr="$chr" ~/mapping/09.filterVCF/GT_AD/field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > ${chr}_119_distance
./maternal_patenal_define.py --son_define on  --i1 ~/mapping/fieldergenomecompare/20200424_201_sample_compare/S140_NZ05/${chr}.homo_hmm_snp_level --i2 ~/mapping/fieldergenomecompare/20200424_201_sample_compare/YM03_S140/${chr}.homo_hmm_snp_level --s1 NZ05 --s2 YM03 -o ${chr}.nd5181_son_belong
./maternal_patenal_define.py --son_define on  --i1 ~/mapping/fieldergenomecompare/20200424_201_sample_compare/YM03_S6554/${chr}.homo_hmm_snp_level --i2 ~/mapping/fieldergenomecompare/20200424_201_sample_compare/YM03_S67/${chr}.homo_hmm_snp_level --s1 S6554 --s2 S67 -o ${chr}.nd3097_son_belong
done
fi

if false;then
for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
for item in {${sample1},${sample2},both,unknown};do

sed '1d'  ${chr}.${son}_son_belong|cut -f1,2,3,6|grep ${item} |sort -n -k 2 > ${item}_${son}_raw.txt 
a=`du -h ${item}_${son}_raw.txt`
b=`echo $a |awk '{print $1}'`
rm ${item}_${son}.txt
if [ "${b}" != 0 ];then
bedtools merge -i ${item}_${son}_raw.txt  -c 4 -o distinct > ${item}_${son}.txt
else
touch ${item}_${son}.txt
fi
done
cat {${sample1},${sample2},both,unknown}_${son}.txt| sort -n -k 2  > ${chr}.${son}_son_norm
./maternal_patenal_define.py --son_norm on  --i1 ${chr}.${son}_son_norm  -o ${chr}.${son}_son_normed
echo ${chr}
./count_CNV.py --split on -i ${chr}.${son}_son_normed -o ${chr}.${son}_son_normed_splited
rm ${chr}.${son}_son_normed
rm ${chr}.${son}_son_norm
cut -f4 ${chr}.${son}_son_normed_splited | sponge ${chr}.${son}_son_normed_splited
cut -f1,2,3,4,5 ${chr}.${son}_son_belong|sort -n -k 2 | paste - ${chr}.${son}_son_normed_splited | sponge ${chr}.${son}_son_normed_splited
done
fi

#sed '1d' ABD_similar_1.csv | xargs echo -n|sed "s/) /\n/g" |sed "s/,//g"|sed "s/ /\t/g"|sed "s/c(//g" >ABD_similar_1_change.csv 
#awk  '{for(i=2;i<=NF;i++) print $1"_"$i}' ABD_similar_1_change.csv|less 
if true;then
for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
cut -f6 ${chr}.nd3097_son_normed_splited|sed  s/"unknown"/"YM03_unknown"/g|sed  s/"both"/"S67_S6554_both"/g | paste ${chr}.nd5181_son_normed_splited -|cut -f1,2,3,6,7|sed  s/"\<unknown\>"/"S140_unknown"/g|sed  s/"\<both\>"/"YM03_NZ05_both"/g > ${chr}_nd3097_nd5181_pedigree
cut -f1,2,3,4 ${chr}_nd3097_nd5181_pedigree| sed '1d' | sponge ${chr}_nd5181_pedigree_final
#awk '{if ($4=="YM03" ){print $1"\t"$2"\t"$3"\t"$5} else {print $1"\t"$2"\t"$3"\t"$4}}' ${chr}_nd3097_nd5181_pedigree > ${chr}_nd5181_pedigree_final
#sed '1d' ${chr}_nd5181_pedigree_final | sponge ${chr}_nd5181_pedigree_final
cut -f1,2,3,5 ${chr}_nd3097_nd5181_pedigree| sed '1d' | sponge ${chr}_nd3097_pedigree_final
done
fi
