#!/usr/bin/env bash

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
echo ($OPTIND-1)
echo $[($OPTIND-1)]

shift $[($OPTIND-1)]
echo "there are $# group files"



for ((i=1; i<=$#; i++)); do
    #echo $i
    printf "%d\t%s\t%s\n" $i "${!i}" ${!i}
    #Group=$(basename "${!i}"|cut -d. -f1)
    #gawk -vG=${Group} '{print $1"\t"G}' "${!i}" 
    #gawk -vG=${Group} '{print $1}' "${!i}"
    #echo $Group
done
