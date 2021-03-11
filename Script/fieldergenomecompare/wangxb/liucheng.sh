#parallel -j procfile sh bisample_compare_phase1_2.sh ::: $(eval cat sample_list.txt) ::: 2020-04-17162314
#parallel -j procfile sh bisample_compare_phase2_inside.sh ::: $(eval cat sample_list.txt)
#parallel -j procfile sh bisample_compare_phase3_inside.sh ::: $(eval cat sample_list.txt)
declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
done < meta_data_all.txt
final_num=10
while read line;do
           #myjob=`jobs|wc -l`
           SAMPLE1=`echo $line | awk -F ',' '{print $1}'`
           SAMPLE2=`echo $line | awk -F ',' '{print $2}'`
           echo -e "${arr_3[${SAMPLE1}]}\t${arr_3[${SAMPLE2}]}\t$final_num"
done < sample_path.txt
