#!/bin/bash
# set -x

for i in $@;do
  gawk -vt=${i} '{print $1"\t"t}' ${i} | sponge ${i}
done

cat $@ > all.clades.txt
gawk 'ARGIND==1{a[NR]=$1;n=length(a)} ARGIND==2{if(FNR==1){printf "Accession\t";for(i=1;i<=n;i++){printf a[i]"\t"};print ""};printf $1"\t";for(i=1;i<=n;i++){printf $2==a[i]?1"\t":0"\t"}print ""}' <(cut -f2 all.clades.txt|uniq) all.clades.txt | sed 's/\t$//' > clade_map_input.txt

#? : conditional expression :awk '{max=($1 > $2) ? $1 : $2;print max}' filename. if $1>$2 max = $1 else max = $2 
