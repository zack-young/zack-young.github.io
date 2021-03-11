#!/usr/bin/env bash
Rscript CNV_heatmap.R -i chr2B.1M_600M_800M_plus_sample.join -M metadata_jmlx.txt -o chr2B.1M_500M_620M_plus_sample
Rscript CNV_heatmap.R -i chr2D.1M_500M_620M_plus_sample.join -M metadata_jmlx.txt -o chr2D.1M_500M_620M_plus_sample
Rscript CNV_heatmap.R -i chr4B.1M_250M_310M_plus_sample.join -M metadata_jmlx.txt -o chr4B.1M_250M_310M_plus_sample


#paste chr2D.1M.norm_lx99 chr2D.1M.norm_jm22 chr2D.1M.join |sed -n '1.501,621p' > chr2D.1M_500M_620M_plus_sample.join
#sed -i "1i`paste chr4B.1M.norm_lx99 chr4B.1M.norm_jm22 chr4B.1M.join |sed -n "1p"`" chr4B.1M_250M_310M_plus_sample.join


