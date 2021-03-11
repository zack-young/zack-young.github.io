#!/usr/bin/env bash
set -x

while getopts 'g:a:f:n:' opt
do
    case "$opt" in
    g) GENE=$OPTARG ;;
    a) LABEL=$OPTARG ;;
    f) FLANK=$OPTARG ;;
    n) NONE=$OPTARG ;;
    esac
done

if [[ -z "${GENE}" ]]; then
    echo "No gene provided."; exit 1
fi
if [[ ! -e "${LABEL}" ]]; then 
    echo "No annotation file provided"; exit 2
fi
if [[ -z "${FLANK}" ]]; then 
    FLANK=0
fi
if [[ ! -e "${NONE}" ]]; then 
    NONE="none-var-gene.txt"
fi

shift $[($OPTIND-1)]
echo "there are $# group files"

tmpfileS=`mktemp -p /data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/`
tmpfileG=`mktemp -p /data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/`
for ((i=1; i<=$#; i++)); do
    #printf "%d\t%s\n" $i "${!i}"
    Group=$(basename "${!i}"|cut -d. -f1)
    gawk -vG=${Group} '{print $1"\t"G}' "${!i}" >> $tmpfileG
    gawk -vG=${Group} '{print $1}' "${!i}" >> $tmpfileS
done

i=$(grep -n $GENE IWGSC-ath-rice.ann | cut -d":" -f1)
CHR=$(gawk -vFS="\t" -vl=$i 'NR==l{print $2}' IWGSC-ath-rice.ann)
START=$(expr $(gawk -vFS="\t" -vl=$i 'NR==l{print $3}' IWGSC-ath-rice.ann) - $FLANK)
END=$(expr $(gawk -vFS="\t" -vl=$i 'NR==l{print $4}' IWGSC-ath-rice.ann) + $FLANK)
ANNO=$(gawk -vFS="\t" -vl=$i 'NR==l{print $6"|"$7"|"$8"|"$9"|"$10}' IWGSC-ath-rice.ann)

tmpfile=`mktemp -p /data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/`
bcftools view /data2/rawdata2/mergeFile/version2/${CHR}.ann.bcf.gz \
    -S $tmpfileS \
    --min-ac=1 --threads 4 \
    -r ${CHR}:${START}-${END} |./VCF2Matrix.py > $tmpfile

#echo "" >> ${GENE}.gene
#echo "# "$GENE >> ${GENE}.gene
#echo "# "$ANNO >> ${GENE}.gene
#bcftools view /data2/rawdata2/mergeFile/version2/${CHR}.ann.bcf.gz \
#    -S $tmpfileS \
#    --min-ac=1 --threads 4 \
#    -r ${CHR}:${START}-${END} \
#    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t\n' \
#    | gawk -vg=${GENE} -F'[\t,]' -vOFS="\t" '{printf $1"\t"$2"\t"$3"\t"$4"\t";for(i=5;i<=NF;i++){if($i~g && $i!~/intergenic_region/){printf $i","}};print ""}' >> ${GENE}.gene


if [ $(wc -l $tmpfile | cut -f1 -d" ") -eq 1 ]
then
    echo $GENE" novar" 
    (echo -e ${GENE}"\t\c"
    for ((i=1; i<=$#; i++)); do
        echo -e "${!i}\t\c"
    done
    echo "" ) >> ${NONE}
else
    Rscript PlotHaplotype.R -i $tmpfile -g $tmpfileG -a $LABEL -n $GENE -r ${CHR}:${START}-${END} 
    # -i infile -g groupfile -a annotationfile -n genename -r region
fi

rm $tmpfile
rm $tmpfileS
rm $tmpfileG
