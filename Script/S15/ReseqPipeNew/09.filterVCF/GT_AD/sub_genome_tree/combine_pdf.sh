for maf in 01 02 03 04 05 1;do
for miss in 05 1 2 3 4 5;do
qpdf --empty --pages {DD,AABBDD,AABB}_setname_maf${maf}_miss${miss}.pdf -- all_maf${maf}_miss${miss}.pdf &
sleep 0.5
done
done
wait
#qpdf can decipher the coded pdf file
