#!/bin/bash
#for VCF in /data2/rawdata2/mergeFile/10.concatVCF/chr??.bcf.gz;do
#for VCF in /data/user/yangzz/mapping/5xresult/*.raw.bcf;do
#for VCF in /data2/rawdata2/mergeFile/other_sample/YangZhengzhao/*snp.filter.final.vcf.gz; do
#for VCF in /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/*bcf.gz; do
#    CHR=`basename $VCF | sed s/.bcf.gz//g`
#    CHR=`basename $VCF | sed s/.snp.filter.final.vcf.gz//g`
#    echo "${VCF}"
#    bcftools view  -v indels -Ov $VCF > ${CHR}.indel.vcf &
#    bcftools view  -v indels -Ov $VCF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > GT_AD/lx99_jm22/${CHR}_lx99_jm22_indel.gt &
    #bcftools view  -S samplename_lx99_jm22.txt -Ov $VCF > ${CHR}.bcf.gz &
#    bcftools view -g '^miss' -i 'FORMAT/GT="AA"' $VCF |bcftools view -e 'FORMAT/DP[0-1]<3 & FORMAT/DP[0-1]>20' > ${CHR}.filter_vcf &
#    bcftools view -S samplename.txt -Ov $VCF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > GT_AD/F1_2_gt/${CHR}_f1_2.gt &
#bcftools view -S samplename.txt -Ov | bcftools query -i'FMT/DP>4 && FMT/GQ>5' -f'%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' > gt$1.txt
#done


for VCF in /data/user/yangzz/mapping/S15/ReseqPipeNew/08.mergeGVCF/chr??.bcf.gz; do
    CHR=`basename $VCF | sed s/.bcf.gz//g`
    echo "${VCF}"
    bcftools view  -v indels -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > GT_AD/${CHR}_indel_akcp.gt &
    bcftools view  -v snps -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > GT_AD/${CHR}_snp_akcp.gt &
    #bcftools view  -v indels -s LX99-1-1 -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > GT_AD/lx99_snp_indel_analysis/${CHR}_indel.gt &
    #bcftools view  -v snps -s LX99-1-1 -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > GT_AD/lx99_snp_indel_analysis/${CHR}_snp.gt &
    #bcftools view  -v indels -s C4-2 -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > GT_AD/jm22_snp_indel_analysis/${CHR}_indel.gt &
    #bcftools view  -v snps -s C4-2 -Ov $VCF| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%GQ\t]\n' > GT_AD/jm22_snp_indel_analysis/${CHR}_snp.gt &
done
