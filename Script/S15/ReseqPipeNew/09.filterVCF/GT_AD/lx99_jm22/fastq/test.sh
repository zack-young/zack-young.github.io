touch test_result
for i in {'CATCGTTTTCTCATCAGATCCA','TCATCGTTTTCTCATCAGATCC','GGCTACGAACTCCGAAGTAG','CAGCCAGATCTGTACCTCAG'};do
     a=`wc -l ~/mapping/08.mergeGVCF/08.localBatch.sh`
     echo "$a $i" >> test_result &
done
