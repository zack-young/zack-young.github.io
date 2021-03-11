#!/bin/bash
set -x

while read g t;do
  bash ./script_own.sh -g $g -f 5000 -t $t -o ${t}.pdf s742.grp ak58.grp
done < genelist.txt
