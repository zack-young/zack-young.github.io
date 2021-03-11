# Python 2
import os
import sys
from multiprocessing import Pool, Process, Manager

if len(sys.argv)<3:
	print("python 1_copy_datafile.py <From> <To>")
	exit()

path = sys.argv[1]

folders= os.listdir(path)

new_data_folder_path = sys.argv[2]

if os.path.exists(new_data_folder_path):
	os.system("rm -rf " + new_data_folder_path)

os.makedirs(new_data_folder_path)

def copy_file(folder):
	if os.path.isdir(path + "/" + folder):
		path_c = path + "/" + folder
		folders_c = os.listdir(path_c)
		if "chr1A.homo_undefined_snp_level" in folders_c:
			os.makedirs(new_data_folder_path + "/" + folder)
			for i in range(0,7):
				for j in ["A", "B", "D"]:
					shell = "ln -s " + path_c + "/chr" + str(i+1) + j + ".1M.density" + " " + new_data_folder_path + "/" + folder + "/"
					os.system(shell)
					shell = "ln -s " + path_c + "/chr" + str(i+1) + j + ".homo_undefined_snp_level" + " " + new_data_folder_path + "/" + folder + "/"
					os.system(shell)
	return 0


proce_pool = Pool(processes = 40)
mulres = []
for folder in folders:
	mulres.append(proce_pool.apply_async(copy_file, args=(folder,)))

proce_pool.close()
proce_pool.join()
				
