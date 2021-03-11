parallel -j procfile sh bisample_compare_phase1_2.sh ::: $(eval cat sample_list.txt) ::: 2020-04-17162314
parallel -j procfile sh bisample_compare_phase2_inside.sh ::: $(eval cat sample_list.txt)
parallel -j procfile sh bisample_compare_phase3_inside.sh ::: $(eval cat sample_list.txt)
