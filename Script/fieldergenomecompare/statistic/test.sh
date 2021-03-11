#numerator=15
#denominator=30
#te=$numerator$denominator
#echo  $numerator $denominator
#echo ${te}
#awk 'BEGIN{printf "%.2f\n",('$numerator'+'$numerator')/'$denominator'}'
DEV_PATH='/data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare'
SAMPLE1='S138'
SAMPLE2='S137'
chr='chr1A'
numerator1=`grep 'low' ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/${chr}.homo_hmm_snp_level|wc -l`
numerator2=`grep 'both_CNV' ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/${chr}.homo_hmm_snp_level|wc -l`
denominator=`cat ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/${chr}.homo_hmm_snp_level | wc -l `
decimal=`awk 'BEGIN{printf "'$DEV_PATH'""\t""%.2f\n",('$numerator1'+'$numerator2')/'$denominator'}'`
echo $decimal

#awk '{print "'$DEV_PATH'""/""'$SAMPLE1'"$1"__""'$SAMPLE2'""_dist"}'  ../metadata_cultivar_final.txt
#for((i=1;i<=1;i++));  
#do   
#echo -en '\n'
#echo 'ed'  
#echo
#done 
