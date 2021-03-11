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

echo "optind is $OPTIND" # -g TAD is two argument
echo "first arg $1,num $#"
shift 1
#shift $[($OPTIND-1)]
echo "$[$OPTIND -1]"
echo "$[($OPTIND-1)]"
echo "first arg $1,num $#"
echo "there are $# group files"

for ((i=1; i<=$#; i++)); do
    #printf "%d\t%s\n" $i "${!i}"
    Group=$(basename "${!i}"|cut -d. -f1)
    echo $Group
#    gawk -vG=${Group} '{print $1"\t"G}' "${!i}" 
#    gawk -vG=${Group} '{print $1}' "${!i}" 
done
