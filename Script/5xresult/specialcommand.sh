cat genotype_count.txt |awk 'ORS=NR%6?"\t":"\n"{print}' |sed 's/|/ /' > genotype_count_line.txt
cat chr2A.200_600M.10DP_GQ.clean |awk '{if(substr($5,1,1)!=substr($5,3,1)){print $0;}}'|head -100
