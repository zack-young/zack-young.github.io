#!/bin/bash
for i in {hete,homo,miss};do
  cat indel/chr*indel_${i}_GQ | awk '{a[$1]+=$2;}END{for(i in a){print i"\t"a[i];}}' -|sort -n > indel/combine_indel_${i}_GQ
  cat indel/chr*indel_${i}_DP | awk '{a[$1]+=$2;}END{for(i in a){print i"\t"a[i];}}' -|sort -n > indel/combine_indel_${i}_DP
done

