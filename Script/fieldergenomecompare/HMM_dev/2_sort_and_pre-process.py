# Python 2
import os
import sys
from multiprocessing import Pool, Process, Manager

if len(sys.argv)<5:
	print("python 2_sort_and_pre-process.py <Folder Path> <Level File Name Extension> <Density File Name Extension> <Output File Name Extension>")
	exit()
  
path = sys.argv[1]

file1_backend = sys.argv[2]
file2_backend = sys.argv[3]

folders= os.listdir(path)

def process_data(folder):
	if os.path.isdir(path + "/" + folder):
		path_c = path + "/" + folder
		for i in range(0,7):
			for j in ["A", "B", "D"]:
				shell = "sort -n -k 2 " + path_c + "/" + "chr" + str(i+1) + j + "." + file1_backend + " > " + path_c + "/" + "chr" + str(i+1) + j + "." + file1_backend + ".sorted"
				#print(shell)
				os.system(shell)
				shell = "sort -n -k 2 " + path_c + "/" + "chr" + str(i+1) + j + "." + file2_backend + " > " + path_c + "/" + "chr" + str(i+1) + j + "." + file2_backend + ".sorted"
				os.system(shell)
				f = open(path_c + "/" + "chr" + str(i+1) + j + "." + file1_backend + ".sorted")
				lines1 = f.readlines()
				f.close()
				f = open(path_c + "/" + "chr" + str(i+1) + j + "." + file2_backend + ".sorted")
				lines2 = f.readlines()
				f.close()
				for k in range(0, len(lines1)):
					line = lines1[k].strip().split("\t")
					if line[3] != "low" and line[3] != "mid" and line[3] != "high":
						line[3] = str(0.001-10)
						line = "\t".join(line)
						line += "\n"
						lines2[k+1] = line
				f = open(path_c + "/" + "chr" + str(i+1) + j + "." + file2_backend + ".sorted." + sys.argv[4], "w")
				f.writelines(lines2)
				f.close()
	return 0


proce_pool = Pool(processes = 40)
mulres = []
for folder in folders:
	mulres.append(proce_pool.apply_async(process_data, args=(folder,)))

proce_pool.close()
proce_pool.join()
