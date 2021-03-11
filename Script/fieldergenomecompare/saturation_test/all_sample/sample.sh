for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
for i in 0.05 0.1 0.25 0.5 0.75 0.85 ;do
samtools view -s ${i} -b /data2/user2/xiexm/pro_ject/resequencing/compare/s30/05.mergeAsplit/${CHR}.1.dedup.bam > ${CHR}.1.dedup${i}_sort.bam &
samtools view -s ${i} -b /data2/user2/xiexm/pro_ject/resequencing/compare/s30/05.mergeAsplit/${CHR}.2.dedup.bam > ${CHR}.2.dedup${i}_sort.bam &
done
wait_all
done
