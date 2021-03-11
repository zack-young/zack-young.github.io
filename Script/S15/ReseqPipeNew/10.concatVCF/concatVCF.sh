#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

source /WORK/app/osenv/ln1/set2.sh

WP="/WORK/pp192/"
CHR=$1

bcftools concat `ls ../09.filterVCF/${CHR}.??.snp.filter.final.vcf.gz` -o ${CHR}.snp.bcf.gz -O b --threads 6 &
bcftools concat `ls ../09.filterVCF/${CHR}.??.indel.filter.final.vcf.gz` -o ${CHR}.indel.bcf.gz -O b --threads 6 &
bcftools concat `ls ../09.filterVCF/${CHR}.??.*.filter.final.vcf.gz` -o ${CHR}.bcf.gz -O b --threads 6 &

#
wait