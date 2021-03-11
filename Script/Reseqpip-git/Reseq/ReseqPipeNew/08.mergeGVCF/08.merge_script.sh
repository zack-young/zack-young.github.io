#!/usr/bin/env bash
set -euxo pipefail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 2400

# there should be changed to samples you want to merge
SMlst=S0,S1

source /WORK/app/osenv/ln1/set2.sh

WP="/WORK/pp192/"

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start merge GVCF"

FSM=`echo ${SMlst}|cut -d, -f1`
FCHR='chr1A.1'

# for CHR in chr1A.1;do
# for CHR in `ls ${WP}/ReseqPJ/${FSM}/07.splitGVCF/*.vcf.gz | cut -d/ -f7 | cut -d. -f1,2 | uniq;do`
for CHR in `zcat ${WP}/ReseqPJ/${FSM}/07.splitGVCF/${FCHR}.aa.g.vcf.gz|grep contig|gawk -F'=|,' '{print $3}'`;do
    echo ${CHR}
    
    # make a new file 
    echo '#!/usr/bin/env sh' > mergeGVCF_${CHR}.sh
    for GVCF in `ls ${WP}/ReseqPJ/${FSM}/07.splitGVCF/${CHR}.??.g.vcf.gz`;do
        ID=`basename $GVCF|cut -d. -f3`
        echo ${ID}
        #
        (echo "java -Xmx3g \\";
        echo "  -Djava.io.tmpdir=${WP}/tmp \\";
        echo "  -jar ${WP}/Install/GenomeAnalysisTK.jar \\";
        echo "  -R ${WP}/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \\";
        echo "  -T GenotypeGVCFs \\";) >> mergeGVCF_${CHR}.sh
        #
        (for SM in `echo ${SMlst} | tr "," " "`; do
            echo "  --variant ${WP}/ReseqPJ/${SM}/07.splitGVCF/${CHR}.${ID}.g.vcf.gz \\";
        done) >> mergeGVCF_${CHR}.sh
        #
        echo "  -o ${CHR}.${ID}.raw.vcf.gz &"  >> mergeGVCF_${CHR}.sh
    done
    echo "wait"  >> mergeGVCF_${CHR}.sh
    yhbatch -n 1 -e ${CHR}.e -o ${CHR}.o ./mergeGVCF_${CHR}.sh
done

echo "["`date +%Y-%m-%d,%H:%M:%S`"] Start merge GVCF"