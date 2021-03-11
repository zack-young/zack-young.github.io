for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do       
#sort  -nk 2,2 ${chr}.homo_hmm_snp_level |cut -f4|sed '1i\LEVEL'  |paste ${chr}.1M.density -|sed s/low/sim/g|sed s/high/diff/g|sed s/mid/diff/g > ${chr}.1M.density_level
parallel -j procfile  ./bisample_compare_hap_length_type_count.sh ::: ${chr} ::: $(eval cat ../sample_path.txt)


#chr="chr1A"
#while read line;do
#lis=`grep ${line} ../statistic/GSR_MCL/sample_path.txt| awk '{print $1"/"chr".homo_undefined_snp_level"}' chr="$chr"| tr '\n' ' '|sed s/' '$//g`
#num=`cat ${lis} | grep -E 'low|both'|cut -f1,2|sort -k1,2|uniq|wc  -l`
#echo $num $line
#done < ../sample_list.txt
done
