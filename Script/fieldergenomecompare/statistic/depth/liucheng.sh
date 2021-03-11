while read line;do
sample=`echo ${line}|awk '{print $1}'`
name=`echo ${line}|awk '{print $3}'`
total_ref=`cut -d '=' -f4 /data/user/yangzz/mapping/08.mergeGVCF/header_length.txt|sed 's/>//g' | awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'`
sample_counts=`grep 'total' /data2/public/ReseqData/workflowFile/${sample}/05.mergeAsplit/chr*.flagstat | cut -d ' ' -f1|cut -d ':' -f2|awk 'BEGIN{sum=0}{sum+=$1}END{print sum*150}'`
depth=`awk 'BEGIN{printf "'$sample'""\t""'$name'""\t""%.2f\n",('$sample_counts'/'$total_ref')}'`
echo $depth
done < ../../metadata_cultivar_199.txt
