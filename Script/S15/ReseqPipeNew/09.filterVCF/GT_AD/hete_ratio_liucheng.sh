#!/bin/bash
#for NUM in {1..3};do            #draw the picture of DP density plot| use combined data no need to specify chromosome
#    Rscript r_script/DP_hist.R -i /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/DP/T${NUM}_DP/combine_data \
#    -p /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/ \
#    -H 6 -W 12 -f pdf -t "Density plot of T${NUM}_deepth" \
#    -o T${NUM}_deepth.pdf -e T${NUM}.txt
        #
#        Rscript r_script/DP_hist.R -i /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/DP/F${NUM}_DP/chr${num}B.DP \
#        -p /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/ \
#        -H 6 -W 12 -f pdf -t "Density plot of F${NUM}_${num}B deepth" \
#        -o F${NUM}_${num}B_deepth.pdf -e F${NUM}_${num}B.txt
#        #
#        Rscript r_script/DP_hist.R -i /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/DP/F${NUM}_DP/chr${num}D.DP \
#        -p /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/ \
#        -H 6 -W 12 -f pdf -t "Density plot of F${NUM}_${num}D deepth" \
#        -o F${NUM}_${num}D_deepth.pdf -e F${NUM}_${num}D.txt
#        #./ratio_filter.py -b 1000000 -m 9.81 -l 9.12 -n S5 -1 5 -i s5_AD/chr1A.gt > filtered_ratio/s5_ratio/chr1A.ratio &
#done
#wait
#R=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/S1_1A.txt | awk '{print $1}'`
#L=cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/s1_1A.txt | awk '{print $2}'
#
#for NUM in {1..15};do                    # use DP threshold to mask DP and CNV then count the hete ratio
#    for num in {1..7};do
#        if [ ! -d "filtered_ratio/s${NUM}_ratio" ]; then
#            mkdir filtered_ratio/s${NUM}_ratio
#        fi
##        Ma=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/S${NUM}.txt | awk '{printf("%.4f\n",$2)}'`
##        La=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/S${NUM}.txt | awk '{printf("%.4f\n",$1)}'`
#        Ma=
#        La=
#        ./ratio_filter.py -b 1000000 -m ${Ma} -l ${La} -n /data2/rawdata2/readDepth/S${NUM} -1 5 -i s${NUM}_AD/chr${num}A.gt > filtered_ratio/s${NUM}_ratio/chr${num}A.ratio &
#        #echo ${Ma}
#        #
#        ./ratio_filter.py -b 1000000 -m ${Ma} -l ${La} -n /data2/rawdata2/readDepth/S${NUM} -1 5 -i s${NUM}_AD/chr${num}B.gt > filtered_ratio/s${NUM}_ratio/chr${num}B.ratio &
#        #
#        ./ratio_filter.py -b 1000000 -m ${Ma} -l ${La} -n /data2/rawdata2/readDepth/S${NUM} -1 5 -i s${NUM}_AD/chr${num}D.gt > filtered_ratio/s${NUM}_ratio/chr${num}D.ratio &
#    done
#    wait
#    cat filtered_ratio/s${NUM}_ratio/*A.ratio > filtered_ratio/s${NUM}_ratio/A_combine_ratio
#    cat filtered_ratio/s${NUM}_ratio/*B.ratio > filtered_ratio/s${NUM}_ratio/B_combine_ratio
#    cat filtered_ratio/s${NUM}_ratio/*D.ratio > filtered_ratio/s${NUM}_ratio/D_combine_ratio
#done
#wait
#
#find filtered_ratio/  -name "A_combine_ratio" -exec cat '{}' \; > filtered_ratio/A_F_S15

#for NUM in {1..3};do    # specify ABD the draw the picture of hete ratio distribution
##    for num in {1..7};do
#    Sample=T${NUM}
#    Rscript r_script/hete_ratio_genome.R -i /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/filtered_ratio/${Sample}_ratio/A_combine_ratio \
#    -p /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/hete_ratio/ \
#    -H 9 -W 16 \
#    -f pdf \
#    -t "FDR Distribution of ${Sample}_A" \
#    -T "Density plot of ${Sample}_A hete_raito" \
#    -o ${Sample}_A_hete.pdf \
#    -O ${Sample}_A_FDR.pdf \
#    -e ${Sample}_A.txt
#    #
#    Rscript r_script/hete_ratio_genome.R -i /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/filtered_ratio/${Sample}_ratio/B_combine_ratio \
#    -p /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/hete_ratio/ \
#    -H 9 -W 16 \
#    -f pdf \
#    -t "FDR Distribution of ${Sample}_B" \
#    -T "Density plot of ${Sample}_B hete_raito" \
#    -o ${Sample}_B_hete.pdf \
#    -O ${Sample}_B_FDR.pdf \
#    -e ${Sample}_B.txt
#    #
#    Rscript r_script/hete_ratio_genome.R -i /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/filtered_ratio/${Sample}_ratio/D_combine_ratio \
#    -p /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/hete_ratio/ \
#    -H 9 -W 16 \
#    -f pdf \
#    -t "FDR Distribution of ${Sample}_D" \
#    -T "Density plot of ${Sample}_D hete_raito" \
#    -o ${Sample}_D_hete.pdf \
#    -O ${Sample}_D_FDR.pdf \
#    -e ${Sample}_D.txt
##    done
#done
##
#for NUM in {1..3};do     # produce the data to visualize chromsome in homo hete DP CNV
#    for num in {1..7};do
#        if [ ! -d "graphy/T${NUM}_chr_graphy" ]; then
#            mkdir graphy/T${NUM}_chr_graphy        
#        fi
#        CNV=/data/user/yangzz/mapping/5xresult/T${NUM}
#        Ma=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$2)}'`
#        La=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$1)}'`
#        Ha=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/hete_ratio/T${NUM}_A.txt | awk '{printf("%.4f\n",$1)}'`
#        ./chrom_graphy.py -b 1000000 -m ${Ma} -l ${La} -n ${CNV} -e ${Ha} -1 5 -i T${NUM}_AD/chr${num}A.gt > graphy/T${NUM}_chr_graphy/chr${num}A.data &
#        #echo ${Ma}
#        #
#        Mb=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$2)}'`
#        Lb=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$1)}'`
#        Hb=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/hete_ratio/T${NUM}_B.txt | awk '{printf("%.4f\n",$1)}'`
#        ./chrom_graphy.py -b 1000000 -m ${Mb} -l ${Lb} -n ${CNV} -e ${Hb} -1 5 -i T${NUM}_AD/chr${num}B.gt > graphy/T${NUM}_chr_graphy/chr${num}B.data &
#        #
#        Md=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$2)}'`
#        Ld=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$1)}'`
#        Hd=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/hete_ratio/T${NUM}_D.txt | awk '{printf("%.4f\n",$1)}'`
#        ./chrom_graphy.py -b 1000000 -m ${Md} -l ${Ld} -n ${CNV} -e ${Hd} -1 5 -i T${NUM}_AD/chr${num}D.gt > graphy/T${NUM}_chr_graphy/chr${num}D.data &
#    done
#    wait
#done
#
for NUM in {1..3};do
    if [ ! -d "MASK/T${NUM}_mask" ]; then
        mkdir MASK/T${NUM}_mask
    fi
    for num in {1..7};do
        CNV=/data/user/yangzz/mapping/5xresult/T${NUM}
        Ma=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$2)}'`
        La=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$1)}'`
        Ha=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/hete_ratio/T${NUM}_A.txt | awk '{printf("%.4f\n",$1)}'`
        ./mask_DP_CNV_hete.py -b 1000000 -m ${Ma} -l ${La} -n ${CNV} -e ${Ha} -1 5 -i T${NUM}_AD/chr${num}A.gt > MASK/T${NUM}_mask/chr${num}A.mask_CNV &
        #echo ${Ma}
        #
        Mb=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$2)}'`
        Lb=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$1)}'`
        Hb=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/hete_ratio/T${NUM}_B.txt | awk '{printf("%.4f\n",$1)}'`
        ./mask_DP_CNV_hete.py -b 1000000 -m ${Mb} -l ${Lb} -n ${CNV} -e ${Hb} -1 5 -i T${NUM}_AD/chr${num}B.gt > MASK/T${NUM}_mask/chr${num}B.mask_CNV &
        #
        Md=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$2)}'`
        Ld=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/DP/T${NUM}.txt | awk '{printf("%.4f\n",$1)}'`
        Hd=`cat /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/hete_ratio/T${NUM}_D.txt | awk '{printf("%.4f\n",$1)}'`
        ./mask_DP_CNV_hete.py -b 1000000 -m ${Md} -l ${Ld} -n ${CNV} -e ${Hd} -1 5 -i T${NUM}_AD/chr${num}D.gt > MASK/T${NUM}_mask/chr${num}D.mask_CNV &
    done
    wait
    echo "T${NUM}"
done
#
#for NUM in {1..3};do
#    ./ptf_maker.py -p /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/graphy/ \
#    -s T${NUM} \
#    -o /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/direct_file/T${NUM}.direct
#done

#for NUM in {1..3};do
#    Rscript r_script/hete_homo_DP_CNV_graphy.R -i /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/direct_file/T${NUM}.direct \
#    -p /data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/plot_box/chromsome/ \
#    -H 10 -W 15.69 \
#    -f pdf \
#    -t "T${NUM}" \
#    -o T${NUM}_chromsome.pdf
#done
