CHRlst=chr1A.1,chr1A.2,chr1B.1,chr1B.2,chr1D.1,chr1D.2,chr2A.1,chr2A.2,chr2B.1,chr2B.2,chr2D.1,chr2D.2,chr3A.1,chr3A.2,chr3B.1,chr3B.2,chr3D.1,chr3D.2,chr4A.1,chr4A.2,chr4B.1,chr4B.2,chr4D.1,chr4D.2,chr5A.1,chr5A.2,chr5B.1,chr5B.2,chr5D.1,chr5D.2,chr6A.1,chr6A.2,chr6B.1,chr6B.2,chr6D.1,chr6D.2,chr7A.1,chr7A.2,chr7B.1,chr7B.2,chr7D.1,chr7D.2,chrUn.1,chrUn.2
for CHR in `echo ${CHRlst} | tr "," " "`; do
samtools view -s 3.01 -b /data2/public/ReseqData/workflowFile/S12/05.mergeAsplit/${CHR}.dedup.bam -o BAM/${CHR}.01_dedup.bam
samtools view -s 3.03 -b /data2/public/ReseqData/workflowFile/S12/05.mergeAsplit/${CHR}.dedup.bam -o BAM/${CHR}.03_dedup.bam
samtools view -s 3.05 -b /data2/public/ReseqData/workflowFile/S12/05.mergeAsplit/${CHR}.dedup.bam -o BAM/${CHR}.05_dedup.bam
samtools view -s 3.07 -b /data2/public/ReseqData/workflowFile/S12/05.mergeAsplit/${CHR}.dedup.bam -o BAM/${CHR}.07_dedup.bam
samtools view -s 3.09 -b /data2/public/ReseqData/workflowFile/S12/05.mergeAsplit/${CHR}.dedup.bam -o BAM/${CHR}.09_dedup.bam
#samtools view -s 3.3 -b /data2/public/ReseqData/workflowFile/S12/05.mergeAsplit/${CHR}.dedup.bam -o ${CHR}.30_dedup.bam
#samtools view -s 3.6 -b /data2/public/ReseqData/workflowFile/S12/05.mergeAsplit/${CHR}.dedup.bam -o ${CHR}.60_dedup.bam
#samtools view -s 3.9 -b /data2/public/ReseqData/workflowFile/S12/05.mergeAsplit/${CHR}.dedup.bam -o ${CHR}.90_dedup.bam
done
