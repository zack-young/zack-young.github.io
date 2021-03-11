for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
vcftools --bcf ~/mapping/08.mergeGVCF/201_final/${CHR}_change_hete2miss.bcf --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01 --keep sample_list_CNL --window-pi 1000000 --out ${CHR}.CNL.1M_nohet &
#vcftools --bcf ~/mapping/08.mergeGVCF/201_final/${CHR}.bcf.gz --keep sample_list --window-pi 1000000 --out ${CHR}_CN &
vcftools --bcf ~/mapping/08.mergeGVCF/201_final/${CHR}_change_hete2miss.bcf --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01 --keep sample_list_CNL --TajimaD 1000000 --out ${CHR}.CNL.1M_nohet &
wait_all
done
