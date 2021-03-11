import re
out = open('chr1A.1.raw.gt','w')
with open('chr1A.1.raw.vcf','r') as mydata:
    out.write('chr'+'\t'+'pos'+'\t'+'ref'+'\t'+'alt'+'\t'+'sample'+'\t'+'sample'+'\t'+'sample'+'\n')
    for line in mydata:
        if re.match('#.*',line[0]):
            continue
        else:
            line = line.strip().split()
            word1 = line[9].strip().split(':')
            word2 = line[10].strip().split(':')
            word3 = line[11].strip().split(':')
            out.write(line[0]+'\t'+line[1]+'\t'+line[3]+'\t'+line[4]+'\t'+word1[0]+'\t'+word2[0]+'\t'+word3[0]+'\t'+'\n')

out.close()
