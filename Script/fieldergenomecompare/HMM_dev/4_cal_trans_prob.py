#python 2
import sys

if len(sys.argv)<2:
	print("python 4_cal_trans_prob.py <Output Folder Path>")
	exit()

all_lines = []

f = open(sys.argv[1] + "/A-strings-raw.txt")
all_lines.extend(f.readlines())
f.close()
f = open(sys.argv[1] + "/B-strings-raw.txt")
all_lines.extend(f.readlines())
f.close()
f = open(sys.argv[1] + "/D-strings-raw.txt")
all_lines.extend(f.readlines())
f.close()

count = 0.0
count_AA = 0.0
count_AB = 0.0
count_AC = 0.0
count_AD = 0.0
count_BA = 0.0
count_BB = 0.0
count_BC = 0.0
count_BD = 0.0
count_CA = 0.0
count_CB = 0.0
count_CC = 0.0
count_CD = 0.0
count_DA = 0.0
count_DB = 0.0
count_DC = 0.0
count_DD = 0.0
count_strange = 0.0

for line in all_lines:
	line = line.strip()
	for i in range(0, len(line)-1):
		if line[i:i+2] == "AA":
			count += 1
			count_AA += 1
		elif line[i:i+2] == "AB":
			count += 1
			count_AB += 1
		elif line[i:i+2] == "AC":
			count += 1
			count_AC += 1
		elif line[i:i+2] == "AD":
			count += 1
			count_AD += 1
		elif line[i:i+2] == "BA":
			count += 1
			count_BA += 1
		elif line[i:i+2] == "BB":
			count += 1
			count_BB += 1
		elif line[i:i+2] == "BC":
			count += 1
			count_BC += 1
		elif line[i:i+2] == "BD":
			count += 1
			count_BD += 1
		elif line[i:i+2] == "CA":
			count += 1
			count_CA += 1
		elif line[i:i+2] == "CB":
			count += 1
			count_CB += 1
		elif line[i:i+2] == "CC":
			count += 1
			count_CC += 1
		elif line[i:i+2] == "CD":
			count += 1
			count_CD += 1
		elif line[i:i+2] == "DA":
			count += 1
			count_DA += 1
		elif line[i:i+2] == "DB":
			count += 1
			count_DB += 1
		elif line[i:i+2] == "DC":
			count += 1
			count_DC += 1
		elif line[i:i+2] == "DD":
			count += 1
			count_DD += 1
		else:
			count_strange += 1

Prob_AA = count_AA / (count_AA+count_AB+count_AC+count_AD)
Prob_AB = count_AB / (count_AA+count_AB+count_AC+count_AD)
Prob_AC = count_AC / (count_AA+count_AB+count_AC+count_AD)
Prob_AD = count_AD / (count_AA+count_AB+count_AC+count_AD)

Prob_BA = count_BA / (count_BA+count_BB+count_BC+count_BD)
Prob_BB = count_BB / (count_BA+count_BB+count_BC+count_BD)
Prob_BC = count_BC / (count_BA+count_BB+count_BC+count_BD)
Prob_BD = count_BD / (count_BA+count_BB+count_BC+count_BD)

Prob_CA = count_CA / (count_CA+count_CB+count_CC+count_CD)
Prob_CB = count_CB / (count_CA+count_CB+count_CC+count_CD)
Prob_CC = count_CC / (count_CA+count_CB+count_CC+count_CD)
Prob_CD = count_CD / (count_CA+count_CB+count_CC+count_CD)

Prob_DA = count_DA / (count_DA+count_DB+count_DC+count_DD)
Prob_DB = count_DB / (count_DA+count_DB+count_DC+count_DD)
Prob_DC = count_DC / (count_DA+count_DB+count_DC+count_DD)
Prob_DD = count_DD / (count_DA+count_DB+count_DC+count_DD)

print "#File: all"
print "#Prob AA: " + str(Prob_AA)
print "#Prob AB: " + str(Prob_AB)
print "#Prob AC: " + str(Prob_AC)
print "#Prob AD: " + str(Prob_AD)
print "#Prob BA: " + str(Prob_BA)
print "#Prob BB: " + str(Prob_BB)
print "#Prob BC: " + str(Prob_BC)
print "#Prob BD: " + str(Prob_BD)
print "#Prob CA: " + str(Prob_CA)
print "#prob CB: " + str(Prob_CB)
print "#prob CC: " + str(Prob_CC)
print "#prob CD: " + str(Prob_CD)
print "#prob DA: " + str(Prob_DA)
print "#prob DB: " + str(Prob_DB)
print "#prob DC: " + str(Prob_DC)
print "#prob DD: " + str(Prob_DD)
print "#Count Strange: " + str(count_strange)

print "\n"

def render_output_prob(curr_p, all_p):
	if curr_p == min(all_p):
		curr_p = 1.0 - (float(str(all_p[0])[:6]) + float(str(all_p[1])[:6]) + float(str(all_p[2])[:6]) + float(str(all_p[3])[:6]) - float(str(min(all_p))[:6]))
		return str(curr_p)
	else:
		return str(curr_p)[:6]

print "AA = " + render_output_prob(Prob_AA, [Prob_AA, Prob_AB, Prob_AC, Prob_AD])
print "AB = " + render_output_prob(Prob_AB, [Prob_AA, Prob_AB, Prob_AC, Prob_AD])
print "AC = " + render_output_prob(Prob_AC, [Prob_AA, Prob_AB, Prob_AC, Prob_AD])
print "AD = " + render_output_prob(Prob_AD, [Prob_AA, Prob_AB, Prob_AC, Prob_AD])
print "BA = " + render_output_prob(Prob_BA, [Prob_BA, Prob_BB, Prob_BC, Prob_BD])
print "BB = " + render_output_prob(Prob_BB, [Prob_BA, Prob_BB, Prob_BC, Prob_BD])
print "BC = " + render_output_prob(Prob_BC, [Prob_BA, Prob_BB, Prob_BC, Prob_BD])
print "BD = " + render_output_prob(Prob_BD, [Prob_BA, Prob_BB, Prob_BC, Prob_BD])
print "CA = " + render_output_prob(Prob_CA, [Prob_CA, Prob_CB, Prob_CC, Prob_CD])
print "CB = " + render_output_prob(Prob_CB, [Prob_CA, Prob_CB, Prob_CC, Prob_CD])
print "CC = " + render_output_prob(Prob_CC, [Prob_CA, Prob_CB, Prob_CC, Prob_CD])
print "CD = " + render_output_prob(Prob_CD, [Prob_CA, Prob_CB, Prob_CC, Prob_CD])
print "DA = " + render_output_prob(Prob_DA, [Prob_DA, Prob_DB, Prob_DC, Prob_DD])
print "DB = " + render_output_prob(Prob_DB, [Prob_DA, Prob_DB, Prob_DC, Prob_DD])
print "DC = " + render_output_prob(Prob_DC, [Prob_DA, Prob_DB, Prob_DC, Prob_DD])
print "DD = " + render_output_prob(Prob_DD, [Prob_DA, Prob_DB, Prob_DC, Prob_DD])
