#!/bin/bash
for BAM in ~/mapping/deepth/*.bam; do
  samtools index $BAM
done
