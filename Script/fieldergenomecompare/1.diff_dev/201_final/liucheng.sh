cat 201_final/header.txt | tr '\t' '\n' | awk '{print $1"\t"NR}' > 201_final/header_t.txt1
gawk 'ARGIND==1{a[$1]=NR} ARGIND==2{if($1 in a){print a[$1]}}' 201_final/header_t.txt 200_2_nearest.txt > 200_2_nearest_line.txt
