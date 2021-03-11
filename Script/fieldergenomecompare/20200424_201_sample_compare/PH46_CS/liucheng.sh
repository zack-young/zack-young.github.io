for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do       
sort  -nk 2,2 ${chr}.homo_hmm_snp_level |cut -f4|sed '1i\LEVEL'  |paste ${chr}.1M.density -|sed s/low/sim/g|sed s/high/diff/g|sed s/mid/diff/g > ${chr}.1M.density_level
done
