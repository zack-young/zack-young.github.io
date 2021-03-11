a=`grep 'total' /data3/user3/publicdata/ReseqData/workflowFile/SRR512/05.mergeAsplit/chr*.2.flagstat | cut -d ' ' -f1|cut -d ':' -f2|awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'`
b=`grep 'total' /data3/user3/publicdata/ReseqData/workflowFile/SRR512/05.mergeAsplit/chr*.1.flagstat | cut -d ' ' -f1|cut -d ':' -f2|awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'`
total_depth=`echo $a $b |awk '{print ($1+$2)*150/14547261565}'`

#for i in  1 3 5 6 7 8 9 11 13;do
i=9
#mkdir sample_depth_${i}
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
sub_depth=`echo $i $total_depth |awk '{printf("%.3f\n",$1 / $2)}'`
echo $sub_depth
samtools view -s ${sub_depth} -b /data3/user3/publicdata/ReseqData/workflowFile/SRR512/05.mergeAsplit/${CHR}.1.dedup.bam > sample_depth_${i}/${CHR}.1.dedup.bam &
samtools view -s ${sub_depth} -b /data3/user3/publicdata/ReseqData/workflowFile/SRR512/05.mergeAsplit/${CHR}.2.dedup.bam > sample_depth_${i}/${CHR}.2.dedup.bam &
wait_all
done
#done
