# Python 2
import os
import sys

if len(sys.argv)<4:
	print("python 3_make_train_data.py <Folder Path> <Sorted Level File Name Extension> <Output Folder Path>")
	exit()
	
	
path = sys.argv[1]

file_backend = sys.argv[2]

folders= os.listdir(path)

for folder in folders:
	if os.path.isdir(path + "/" + folder):
		path_c = path + "/" + folder
		folders_c= os.listdir(path_c)
		for folder_c in folders_c:
			if folder_c[4:] == "A."+file_backend:
				f = open(path_c + "/" + folder_c)
				lines = f.readlines()
				f.close()
				string_lv = ""
				for line in lines:
					line_ele = line.strip().split("\t")
					if line_ele[3] == "low":
						string_lv += 'C'
					elif  line_ele[3] == "mid":
						string_lv += 'B'
					elif  line_ele[3] == "high":
						string_lv += 'A'
					else:
						string_lv += 'D'
				f = open(sys.argv[3] + "/A-strings-raw.txt", 'a')
				print >> f, string_lv
				f.close()
			if folder_c[4:] == "B."+file_backend:
				f = open(path_c + "/" + folder_c)
				lines = f.readlines()
				f.close()
				string_lv = ""
				for line in lines:
					line_ele = line.strip().split("\t")
					if line_ele[3] == "low":
						string_lv += 'C'
					elif  line_ele[3] == "mid":
						string_lv += 'B'
					elif  line_ele[3] == "high":
						string_lv += 'A'
					else:
						string_lv += 'D'
				f = open(sys.argv[3] + "/B-strings-raw.txt", 'a')
				print >> f, string_lv
				f.close()
			if folder_c[4:] == "D."+file_backend:
				f = open(path_c + "/" + folder_c)
				lines = f.readlines()
				f.close()
				string_lv = ""
				for line in lines:
					line_ele = line.strip().split("\t")
					if line_ele[3] == "low":
						string_lv += 'C'
					elif  line_ele[3] == "mid":
						string_lv += 'B'
					elif  line_ele[3] == "high":
						string_lv += 'A'
					else:
						string_lv += 'D'
				f = open(sys.argv[3] + "/D-strings-raw.txt", 'a')
				print >> f, string_lv
				f.close()
