f = open("/data/user/shinyug/python/HMMinuse/data/raw/combine5M", "r")
lines = f.readlines()
f.close()
dic = {"lx987":"1", "987":"1", "3097":"2", "hete":"0", "NA":"3"}
par_lis = []
for line in lines:
   col = line.strip().split("\t")
   if col[0].find(".1") != -1 :
       par_lis.append("\t".join([col[0], dic[col[4]]]))
lis = []
curr_chr = ""
curr_ind = -1
for line in par_lis:
    col = line.split("\t")
    if col[0] == curr_chr:
        lis[curr_ind].append(col[1])
    else:
        curr_chr = col[0]
        curr_ind += 1
        lis.append([])
        lis[curr_ind].append(col[1])
len(lis)
for j in range(0,7):
    for i in range(0,3):
        lett = ["A", "B", "D"]
        filename = "/data/user/shinyug/python/HMMinuse/data/5M_parents/7x3/chr" + str(j+1) + lett[i]
        f = open(filename, "w")
        for k in lis[j*3+i]:
            f.write(k)
        f.close()

for i in range(0,3):
    lett = ["A", "B", "D"]
    filename = "/data/user/shinyug/python/HMMinuse/data/5M_parents/1x3/chrx" + lett[i]
    f = open(filename, "w")
    for j in range(0,7):
        for k in lis[j*3 + i]:
            f.write(k)
    f.close()
