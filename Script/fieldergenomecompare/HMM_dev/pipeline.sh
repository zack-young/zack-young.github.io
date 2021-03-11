###########################################
# Note: All paths should NOT end with "/" #
###########################################

# Copy data
# python 1_copy_datafile.py <From> <To>
python 1_copy_datafile.py /data/user/yangzz/mapping/fieldergenomecompare/20200424_201_sample_compare ./201_sample_compare_HMMdata

# Sort and convert
# python 2_sort_and_pre-process.py <Folder Path> <Level File Name Extension> <Density File Name Extension> <Output File Name Extension>
python 2_sort_and_pre-process.py ./201_sample_compare_HMMdata homo_undefined_snp_level 1M.density processed

# (Optional) (To calculate transition matrix)
# Make train data
# python 3_make_train_data.py <Folder Path> <Sorted Level File Name Extension> <Output Folder Path>
python 3_make_train_data.py ./201_sample_compare_HMMdata homo_undefined_snp_level.sorted .

# (Optional) 
# Calculate transition matrix
# python 4_cal_trans_prob.py <Output Folder Path>
python 4_cal_trans_prob.py . > 4_cal_trans_prob.log

# (If the transition matrix was recalculated, the result need to be copied from "4_cal_trans_prob.log" to "5_remake_levelfile_multiproc_ownfunc.py" manually)
# Use GaussianHMM to remake level by converted data file made in step.2
# python3 5_remake_levelfile_multiproc_ownfunc.py <Data Folder Path> <Converted Data File Name Extension> <Output File Name Extension>
/data3/user3/wangwx/bin/Python3.8/bin/python3 5_remake_levelfile_multiproc_ownfunc.py ./201_sample_compare_HMMdata 1M.density.sorted.processed homo_undefined_snp_level.sorted.v1 > 5_remake_levelfile_multiproc_ownfunc.log

# (Optional) 
# Check if any nonCNVs were wrongly regarded as CNVs
# python 6_check_level_2_CNV.py <Remake Log Path>
python 6_check_level_2_CNV.py 5_remake_levelfile_multiproc_ownfunc.log