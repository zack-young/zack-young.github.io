#!/usr/bin/env bash
set -x

while getopts 'c:s:e:a:f:n:p:' opt
do
    case "$opt" in
    c) CHR=$OPTARG ;;
    s) START=$OPTARG ;;
    e) END=$OPTARG ;;
    a) LABEL=$OPTARG ;;
    f) FLANK=$OPTARG ;;
    n) NONE=$OPTARG ;;
    p) path=$OPTARG ;;
    esac
done

if [[ -z "${CHR}" ]]; then
    echo "No chromosome provided."; exit 1
fi
if [[ -z "${START}" ]]; then
    echo "No start_position provided."; exit 1
fi
if [[ -z "${END}" ]]; then
    echo "No end_position provided."; exit 1
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

#ANNO=$(gawk -vFS="\t" -vl=$i 'NR==l{print $6"|"$7"|"$8"|"$9"|"$10}' IWGSC-ath-rice.ann)
if false;then
tmpfile=`mktemp -p /data/user/yangzz/mapping/09.filterVCF/GT_AD/snp_heatmap/`
#bcftools view /data2/rawdata2/mergeFile/version2/${CHR}.ann.bcf.gz \
#bcftools view ${path}/${CHR}.ann.bcf.gz \
bcftools view /data2/rawdata2/shinyKiloProj/KiloProj.ann.bcf.gz \
    --min-af=0.05 \
    -m 2 \
    -M 2 \
    -v snps \
    -S $tmpfileS \
    --min-ac=1 --threads 4 \
    -r ${CHR}:${START}-${END} \
    | ./VCF2Matrix_nofilter.py -c  "-1,0,0.5,1"> $tmpfile
#    | ./VCF2Matrix.py -c  "-1,0,0.5,1"> $tmpfile

#echo "" >> ${GENE}.gene
#echo "# "$GENE >> ${GENE}.gene
#echo "# "$ANNO >> ${GENE}.gene
#bcftools view /data2/rawdata2/mergeFile/version2/${CHR}.ann.bcf.gz \
#    -S $tmpfileS \
#    --min-ac=1 --threads 4 \
#    -r ${CHR}:${START}-${END} \
#    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t\n' \
#    | gawk -vg=${GENE} -F'[\t,]' -vOFS="\t" '{printf $1"\t"$2"\t"$3"\t"$4"\t";for(i=5;i<=NF;i++){if($i~g && $i!~/intergenic_region/){printf $i","}};print ""}' >> ${GENE}.gene


#Rscript plot_hclust.R -i $tmpfile -g $tmpfileG -a $LABEL -n "n" -r ${CHR}:${START}-${END} -o ${CHR}_${START}_${END}
#Rscript plot_tree_heatmap.R -i $tmpfile -g $tmpfileG -a $LABEL -n "n" -r ${CHR}:${START}-${END} -o ${CHR}_${START}_${END}
    # -i infile -g groupfile -a annotationfile -n genename -r region #"$ANNO"

echo "tmpfileS" $tmpfile
echo 'tmpfileG' $tmpfileG
echo 'tmpfile' $tmpfile
fi
#rm $tmpfile
#rm $tmpfileS
#rm $tmpfileG
