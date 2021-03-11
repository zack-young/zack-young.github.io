for i in `seq 0 9`;do
#echo $(( 18100000 + $i*100000 ))
sh script_region.sh -a ../../metadata_cultivar_final_headed.txt  -c chr4B -s $(( 582000000 + $i*100000 )) -e $(( 582100000 + $i*100000 )) -p ~/mapping/08.mergeGVCF/201_final   use_samplelist
done
