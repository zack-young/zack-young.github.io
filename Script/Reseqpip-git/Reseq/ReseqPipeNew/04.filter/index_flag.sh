#!/usr/bin/env bash
set -euxo pipefail
for BAM in $1/??.filter.bam; do
  echo $BAM
  ID=`basename ${BAM} | sed s/.filter.bam//g`
  samtools index ${BAM} > $1/${ID}.filter.bam.bai &
  samtools flagstat ${BAM} > $1/${ID}.flagstat &
#   rm ${BAM})&
done

