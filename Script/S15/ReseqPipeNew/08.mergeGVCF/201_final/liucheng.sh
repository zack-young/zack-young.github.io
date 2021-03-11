for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#bcftools view -v snps --min-ac=1 -M2 -m2  ${CHR}.bcf.gz  |bcftools query -f "[%GT\t]\n"|wc -l >> filteredsnp_summary &
#bcftools view ${CHR}.bcf.gz -v snps|sed "s/0\/1/\.\/\./g" |bcftools view - -Ob -o ${CHR}_change_hete2miss.bcf &
#CHR="chr5A"
(
(bcftools view -h ${CHR}.bcf.gz | grep -v "##contig" | sed "/^##reference/ r header_split.txt";
bcftools view --threads 4 -H -v snps -q 0.01:minor --min-ac=1 -M2 -m2 -r ${CHR}:1-400000000 ${CHR}.bcf.gz|sed "s/0\/1/\.\/\./g" |gawk -vOFS="\t" '{$1=$1".1";print}')|bcftools view --threads 4 -Oz -o ${CHR}_change_hete2miss_1.vcf.gz
(bcftools view -h ${CHR}.bcf.gz | grep -v "##contig" | sed "/^##reference/ r header_split.txt";
bcftools view --threads 4 -H -v snps -q 0.01:minor --min-ac=1 -M2 -m2 -r ${CHR}:400000000-900000000 ${CHR}.bcf.gz|sed "s/0\/1/\.\/\./g" |gawk -vOFS="\t" '{$1=$1".2";$2=$2-400000000;print}')|bcftools view --threads 4 -Oz -o ${CHR}_change_hete2miss_2.vcf.gz
bcftools index -t ${CHR}_change_hete2miss_1.vcf.gz
bcftools index -t ${CHR}_change_hete2miss_2.vcf.gz) &
wait_all
done
wait

CHR="chr5A"
(bcftools view -h ${CHR}.ann.bcf.gz | grep -v "##contig" | sed "/^##reference/ r header_split.txt";
bcftools view --threads 4 -H -v snps -q 0.01:minor --min-ac=1 -M2 -m2 -r ${CHR}:1-400000000 ${CHR}.ann.bcf.gz|sed "s/0\/1/\.\/\./g" |gawk -vOFS="\t" '{$1=$1".1";print}')|bcftools view --threads 4 -Oz -o ${CHR}_change_hete2miss_1.vcf.gz
(bcftools view -h ${CHR}.ann.bcf.gz | grep -v "##contig" | sed "/^##reference/ r header_split.txt";
bcftools view --threads 4 -H -v snps -q 0.01:minor --min-ac=1 -M2 -m2 -r ${CHR}:400000000-900000000 ${CHR}.ann.bcf.gz|sed "s/0\/1/\.\/\./g" |gawk -vOFS="\t" '{$1=$1".2";$2=$2-400000000;print}')|bcftools view --threads 4 -Oz -o ${CHR}_change_hete2miss_2.vcf.gz
bcftools index -t ${CHR}_change_hete2miss_1.vcf.gz
bcftools index -t ${CHR}_change_hete2miss_2.vcf.gz
