#for i in {'CATCGTTTTCTCATCAGATCCA','TCATCGTTTTCTCATCAGATCC','GGCTACGAACTCCGAAGTAG','CAGCCAGATCTGTACCTCAG'};do
for i in {'CATCGTTTTCTCATCAGATCCA','TGGATCTGATGAGAAAACGATG'};do
    touch "${i}.txt"
    echo ${i} >> ${i}.txt
    zcat LX99-1-1_H3TFVDMXX_L1_1.clean.fq.gz | grep ${i} >> ${1}.txt &
    zcat LX99-1-1_H3TFVDMXX_L1_2.clean.fq.gz | grep ${i} >> ${1}.txt &
done
#wait 
#cat *txt > result/total.txt

