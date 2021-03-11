#!/bin/bash
set -euxo pipefail

G1=PPD1-2A
G2=PPD1-2D
WD=/data2/rawdata2/HapNet/3pathways/cycle/

comm -12 <(cut -f 1 ${WD}/${G1}/hap1|sort) <(cut -f 1 ${WD}/${G2}/hap1|sort) > hap3
comm -12 <(cut -f 1 ${WD}/${G1}/hap1|sort) <(cut -f 1 ${WD}/${G2}/hap2|sort) > hap2
comm -12 <(cut -f 1 ${WD}/${G1}/hap2|sort) <(cut -f 1 ${WD}/${G2}/hap1|sort) > hap1
comm -12 <(cut -f 1 ${WD}/${G1}/hap2|sort) <(cut -f 1 ${WD}/${G2}/hap2|sort) > hap0

bash /data2/rawdata2/altitude/lib/turn_file_to_mat.sh hap0 hap1 hap2 hap3
