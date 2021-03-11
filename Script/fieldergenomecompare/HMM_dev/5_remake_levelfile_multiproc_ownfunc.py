#python3
import os
import sys
import math
import numpy
from scipy import stats
from hmmlearn import hmm
from scipy.stats import norm
from multiprocessing import Pool, Process, Manager

if len(sys.argv)<4:
	print("python3 5_remake_levelfile_multiproc_ownfunc.py <Data Folder Path> <Converted Data File Name Extension> <Output File Name Extension>")
	exit()
	
path = sys.argv[1]

def proc_func(folder,states,start_arr,trans_mat,means,covars,weight,cdfne,ofne):
	file1_backend = "homo_undefined_snp_level.sorted"
	file2_backend = cdfne

	new_file_backend = ofne

	log_s = ""

	if os.path.isdir(path + "/" + folder):
		path_c = path + "/" + folder

		for i in range(0,7):
			for j in ["A", "B", "D"]:
				f = open(path_c + "/" + "chr" + str(i+1) + j + "." + file1_backend)
				lines1 = f.readlines()
				f.close()

				f = open(path_c + "/" + "chr" + str(i+1) + j + "." + file2_backend)
				lines2 = f.readlines()
				f.close()
				lines2 = lines2[1:]

				raw_seq = []
				for k in range(len(lines2)):
					line = lines2[k].strip().split("\t")
					raw_seq.append(math.log10(float(line[3])+10))

				# viterbi
				mat = [[0 for m in range(len(raw_seq))] for n in range(4)]
				matTB = [[0 for m in range(len(raw_seq))] for n in range(4)]
				# Fill in first column
				for m in range(0, len(states)):
					pdf_val = norm.pdf(raw_seq[0], loc=means[m], scale=covars[m])
					if pdf_val == 0.0:
						eq = -64
					else:
						eq = math.log10(pdf_val)
					mat[m][0] = start_arr[m] + eq + weight[m]
				# Fill in the rest of mat, and choose elements of matTB
				for n in range(1, len(raw_seq)):
					for m in range(0, len(states)):
						pdf_val = norm.pdf(raw_seq[n], loc=means[m], scale=covars[m])
						if pdf_val == 0.0:
							eq = -64
						else:
							eq = math.log10(pdf_val)
						mx, mxi = mat[0][n-1] + trans_mat[0][m] + eq + weight[m], 0
						for m_former in range(1, len(states)):
							pr = mat[m_former][n-1] + trans_mat[m_former][m] + eq + weight[m]
							if pr > mx:
								mx, mxi = pr, m_former
						mat[m][n], matTB[m][n] = mx, mxi
				# Find the final state which has the maximal probability
				omx, omxi = mat[0][len(raw_seq)-1], 0
				for m in range(1, len(states)):
					if mat[m][len(raw_seq)-1] > omx:
						omx, omxi = mat[m][len(raw_seq)-1], m
				# Trace for the state sequence
				m, p = omxi, [omxi]
				for n in range(len(raw_seq)-1, 0, -1):
					m = matTB[m][n]
					p.append(m)
				
				logprob = omx
				log_s += (path_c + "/" + "chr" + str(i+1) + j) + "\n"
				log_s += str(logprob) + "\n"
				count_lev_to_CNV = 0
				for k in range(len(lines1)):
					state = states[p[len(lines1)-k-1]]
					line = lines2[k].strip().split("\t")
					if float(line[3]) > 0:
						if state == "CNV":
							count_lev_to_CNV += 1
						else:
							line[3] = state
						line = "\t".join(line)
						line += "\n"
						lines1[k] = line
				log_s += str(count_lev_to_CNV) + "\n"
				f = open(path_c + "/" + "chr" + str(i+1) + j + "." + new_file_backend, "w")
				f.writelines(lines1)
				f.close()
	return log_s


#model.transmat_ = model.transmat_ / modelD.transmat_.sum(axis=1)[:, numpy.newaxis]
#model.emissionprob_ = model.emissionprob_ / modelD.emissionprob_.sum(axis=1)[:, numpy.newaxis]

AA = 0.8127
AB = 0.1399
AC = 0.0008
AD = 0.0466
BA = 0.0465
BB = 0.9345
BC = 0.0038
BD = 0.0152
CA = 0.0056
CB = 0.0884
CC = 0.8856
CD = 0.0204
DA = 0.2023
DB = 0.1915
DC = 0.0098
DD = 0.5964

AA = math.log10(AA)
AB = math.log10(AB)
AC = math.log10(AC)
AD = math.log10(AD)
BA = math.log10(BA)
BB = math.log10(BB)
BC = math.log10(BC)
BD = math.log10(BD)
CA = math.log10(CA)
CB = math.log10(CB)
CC = math.log10(CC)
CD = math.log10(CD)
DA = math.log10(DA)
DB = math.log10(DB)
DC = math.log10(DC)
DD = math.log10(DD)

# ["A", "B", "C", "D"]
states = ["high", "mid", "low", "CNV"]
n_states = 4
n_obs = 4

start_arr_sing = math.log10(0.25)
start_arr = [start_arr_sing, start_arr_sing, start_arr_sing, start_arr_sing]
trans_mat = [[AA,AB,AC,AD],[BA,BB,BC,BD],[CA,CB,CC,CD],[DA,DB,DC,DD]]
means = [3.5, 2.1, 1.1, -3]
covars = [0.21149496, 0.31481414, 0.09464676, 0.04473012]
weight = [math.log10(0.1167319), math.log10(0.6282251), math.log10(0.2550430), math.log10(1)]


folders= os.listdir(path)

proce_pool = Pool(processes = 40)
mulres = []
for folder in folders:
	mulres.append(proce_pool.apply_async(proc_func, args=(folder,states,start_arr,trans_mat,means,covars,weight,sys.argv[2],sys.argv[3])))

proce_pool.close()
proce_pool.join()

for res in mulres:
	print(res.get()[:-1])
