#!/usr/bin/env sh

parallel -j 5 Rscript CNV_heat_map_bySample.R -i ../{}/{}.1M.norm.chr -k ../CS_karyotype.txt -p ../pdf/{}_DPbySample.pdf :::: ../samples.txt
