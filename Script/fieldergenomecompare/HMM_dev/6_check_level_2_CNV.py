import sys
if len(sys.argv)<2:
	print("python 6_check_level_2_CNV.py <Remake Log Path>")
	exit()
	
f = open(sys.argv[1])
lines = f.readlines()
f.close()

for i in range(0,len(lines)/3):
    line = lines[i*3+2]
    line = int(line.strip())
    if line>0:
        print(i*3)
        print(line[i*3].strip())
        print(line)
