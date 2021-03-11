#!/usr/bin/env bash
set -xeuo pipefail

while getopts 'g:f:t:o:' opt
do
    case "$opt" in
    g) GENE=$OPTARG ;;
    f) FLANK=$OPTARG ;;
    t) TITLE=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    esac
done

if [ -z "${GENE+x}" ]; then echo "No gene provided."; exit 1; fi
if [ -z "${FLANK+x}" ]; then FLANK=0; fi
if [ -z "${TITLE+x}" ]; then TITLE="Snp Frequency of "${GENE}; fi
if [ -z "${OUTPUT+x}" ]; then OUTPUT=${GENE}.pdf; fi

shift $[($OPTIND-1)]
echo "there are $# group files"

tmpfileS='/data/user/yangzz/tmp/tmpfileS'
tmpfileG='/data/user/yangzz/tmp/tmpfileG'
tmp_group_name='/data/user/yangzz/tmp/tmp_group_name'
for ((i=1; i<=$#; i++)); do
    #printf "%d\t%s\n" $i "${!i}"
    Group=$(basename "${!i}"|cut -d. -f1)
    gawk -vG=${Group} '{print $1"\t"G}' "${!i}" >> $tmpfileG
    gawk -vG=${Group} '{print $1}' "${!i}" >> $tmpfileS
    echo $Group >> $tmp_group_name
done

i=$(grep -n $GENE /data2/rawdata2/database/genome_func_anno/IWGSC-ath-rice.ann | cut -d":" -f1)
CHR=$(gawk -vOFS="\t" -vl=$i 'NR==l{print $2}' /data2/rawdata2/database/genome_func_anno/IWGSC-ath-rice.ann)
START=$(expr $(gawk -vFS="\t" -vl=$i 'NR==l{print $3}' /data2/rawdata2/database/genome_func_anno/IWGSC-ath-rice.ann) - $FLANK)
END=$(expr $(gawk -vFS="\t" -vl=$i 'NR==l{print $4}' /data2/rawdata2/database/genome_func_anno/IWGSC-ath-rice.ann) + $FLANK)
ANNO=$(gawk -vFS="\t" -vl=$i 'NR==l{print $6"|"$7"|"$8"|"$9"|"$10}' /data2/rawdata2/database/genome_func_anno/IWGSC-ath-rice.ann)

tmpbcf='/data/user/yangzz/tmp/tmpbcf.bcf.gz'
tmpfile='/data/user/yangzz/tmp/tmpfile'
mkdir /data/user/yangzz/tmp/tmpsort
tmpsort='/data/user/yangzz/tmp/tmpsort'
bcftools view -M 2 -m 2 -c 2 -r ${CHR}:${START}-${END} -S $tmpfileS --threads 4 /data2/rawdata2/mergeFile/SnpEffAnno/${CHR}.ann.bcf.gz -Ou | bcftools sort -T $tmpsort -Ob -o $tmpbcf
bcftools index $tmpbcf
bcftools view -r ${CHR}:${START}-${END} -S <(cat $1 |cut -f 1) --threads 4 $tmpbcf |bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ANN\t%AC\t%AN\n" > $tmpfile
#for ((i=2; i<=$#; i++)); do
#  paste <(cat $tmpfile) <(bcftools view -r ${CHR}:${START}-${END} -S <(cat ${!i}|cut -f 1) --threads 4 $tmpbcf |bcftools query -f "%AC\t%AN\n") | sponge $tmpfile
#done
gawk -F"\t" -vOFS="\t" -vg=$GENE '{print $0"\t"g}' $tmpfile | sponge $tmpfile

n_group=$(wc -l < $tmp_group_name)
infile='/data/user/yangzz/tmp/infile'
gawk -F"\t" -vOFS="\t" 'ARGIND==1{a=a"\t"$1} END{print "CHR\tPOS\tTYPE"a}' $tmp_group_name >> $infile
gawk -F"\t" -vOFS="\t" '{split($5,anns,",");l=length(anns);if(l==1){ann=anns[1]}else{for(i=1;i<=l;i++){if(anns[i]~$10){ann=anns[i]}}};split(ann,a,"|");if($7==0){afa=0}else{afa=$6/$7};print $1,$2,a[2],afa;delete anns;delete a;ann=""}' $tmpfile >> $infile
# three groups
### gawk -vOFS="\t" '{split($5,anns,",");l=length(anns);if(l==1){ann=anns[1]}else{for(i==1;i<=l;i++){if(anns[i]~$10){ann=anns[i]}}};split(ann,a,"|");print $1,$2,a[2],$6/$7,$8/$9,$10/$11}' tmp.res2 > tmp.res3
gawk -F"\t" -vOFS="\t" '{if(NR>1){if($3~/missense_variant/){$3=1}else if($3~/synonymous_variant/){$3=2}else if($3~/frameshift_variant/){$3=3}else if($3~/stop_/){$3=4}else if($3~/splice_region_variant/){$3=5}else{$3=6}};print}' $infile | sponge $infile 

Rscript lollipop.R -i ${infile} -o ${OUTPUT} -c ${CHR} -l ${START} -r ${END} -t "${TITLE}" -p "${ANNO}"

#rm $tmpfile
#rm $tmpfileS
#rm $tmpfileG
#rm $tmp_group_name
#rm $tmpbcf
#rm $infile
