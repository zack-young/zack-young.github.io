#!/bin/env bash

set -euxo pipefail
for SM in *.vcf; do
   FD=`basename ${SM}|cut -d . -f1`
   ~/tool/annovar/convert2annovar.pl -format vcf4 ${SM}>${FD}.avinput
   perl ~/tool/annovar/table_annovar.pl ./${FD}.avinput ~/genome/IW --outfile ./${FD} --buildver IW --protocol refGene --remove --csvout --operation g
   rm -rf ${FD}.avinput
done
