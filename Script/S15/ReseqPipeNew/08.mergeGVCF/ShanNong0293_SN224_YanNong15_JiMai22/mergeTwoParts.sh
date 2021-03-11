#!/usr/bin/env bash
set -euxo pipefail

CHR=$1
(bcftools view -h ${CHR}.1.sort.bcf.gz | grep -v "##contig" | sed "/^##reference/ r header_length.txt";
bcftools view --threads 4 -H ${CHR}.1.sort.bcf.gz|gawk -vOFS="\t" '{$1=gensub(/(chr..)../,"\\1","g",$1);print}'; 
bcftools view --threads 4 -H ${CHR}.2.sort.bcf.gz|gawk -vOFS="\t" '{$1=gensub(/(chr..)../,"\\1","g",$1);$2=$2+400000000;print}' )|\
bcftools view --threads 4 -Ob -o ${CHR}.bcf.gz

bcftools index ${CHR}.bcf.gz


#Not all versions of sed understand \t. Just insert a literal tab instead (press Ctrl-V then Tab).
# when using sed to replace chr name, it seems that sometimes will fail(has sth to do with TAB)
#(cat header; bcftools view --threads 4 -H chr1A.1.snp.vcf.gz|gawk -vOFS="\t" '{$1=gensub(/(chr..)../,"\\1","g",$1);print}'; bcftools view -H --threads 4 chr1A.2.snp.vcf.gz|gawk -vOFS="\t" '{$1=gensub(/(chr..)../,"\\1","g",$1);$2=$2+400000000;print}' )|bcftools view --threads 4 -Ob -o chr1A.snp.bcf.gz

# substitute header
#(bcftools view -h some_variants.bcf 
#    | grep -v "##contig" 
#    | sed "/^##reference/ r header_length.txt") 
#    | bcftools view -Oz -o contig_1.vcf.gz

# non-split-header save as header_length.txt

(bcftools view -h ${CHR}.1.indel.filter.final.vcf.gz | grep -v "##contig" | sed "/^##reference/ r header_length.txt";
bcftools view --threads 4 -H ${CHR}.1.indel.filter.final.vcf.gz|gawk -vOFS="\t" '{$1=gensub(/(chr..)../,"\\1","g",$1);print}'; 
bcftools view --threads 4 -H ${CHR}.2.indel.filter.final.vcf.gz|gawk -vOFS="\t" '{$1=gensub(/(chr..)../,"\\1","g",$1);$2=$2+400000000;print}' )|\
bcftools view --threads 4 -Ob -o ${CHR}.indel.bcf.gz

bcftools index ${CHR}.indel.bcf.gz

(bcftools view -h ${CHR}.1.snp.filter.final.vcf.gz | grep -v "##contig" | sed "/^##reference/ r header_length.txt";
bcftools view --threads 4 -H ${CHR}.1.snp.filter.final.vcf.gz|gawk -vOFS="\t" '{$1=gensub(/(chr..)../,"\\1","g",$1);print}'; 
bcftools view --threads 4 -H ${CHR}.2.snp.filter.final.vcf.gz|gawk -vOFS="\t" '{$1=gensub(/(chr..)../,"\\1","g",$1);$2=$2+400000000;print}' )|\
bcftools view --threads 4 -Ob -o ${CHR}.snp.bcf.gz

bcftools index ${CHR}.snp.bcf.gz

