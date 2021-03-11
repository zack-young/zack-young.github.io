for num in {1..7}; do
    for a in {A,B,D}; do
        for i in {high,mid}; do 
            #awk '{print $1"\t"$2"\t"$2}' chr${num}${a}_${i}.ann > chr${num}${a}_${i}.bed
            java -jar  /home/wangzh/bin/snpEff_latest_core/snpEff/snpEff.jar -i bed IWGSCv1p1 chr${num}${a}_${i}.bed -no-downstream -no-upstream > chr${num}${a}_${i}_exon_intron.bed
        done
    done
done
